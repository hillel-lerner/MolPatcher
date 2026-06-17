"""
Parses CHARMM forcefield files (.prm, .str), extracts specific topological parameters
for newly formed molecular junctions, converts them to GROMACS units, and writes
them to an included .itp file.
"""

import os


class ForceField:
    """
    Creates a searchable database of CHARMM parameters for dynamic GROMACS topology generation.
    """

    def __init__(self, ff_files: list) -> None:
        """
        Accepts a list of file paths for multiple CHARMM parameter files to be loaded into a searchable database.

        :param ff_files: A list of absolute or relative paths to CHARMM forcefield files.
        :type ff_files: list
        :return: None
        """

        self.ff_files = ff_files

        self.parameters = {
            "atomtypes": {},
            "bonds": {},
            "angles": {},
            "dihedrals_proper": {},
            "dihedrals_improper": {},
        }

        self.junction_output = {k: set() for k in self.parameters}

    @staticmethod
    def find_ff_type(ff_file):
        """
        Determines the type of CHARMM forcefield file based on its extension.

        :param ff_file: The path or name of the forcefield file.
        :type ff_file: str
        :return: The file extension type ('prm', 'str', or 'itp'), or None if unrecognized.
        :rtype: str or None
        """
        base_name = os.path.basename(ff_file)
        if ".prm" in base_name:
            return "prm"
        elif ".str" in base_name:
            return "str"
        elif ".itp" in base_name:
            return "itp"
        return None

    def read_database(self):
        """
        Loops through all provided files and builds the searchable parameter dictionaries.

        :return: None. Populates the internal self.parameters dictionary.
        """

        charmm_mapping = {
            "ATOMS": "atomtypes",
            "BONDS": "bonds",
            "ANGLES": "angles",
            "DIHEDRALS": "dihedrals_proper",
            "IMPROPER": "dihedrals_improper",
        }

        for file in self.ff_files:
            filetype = self.find_ff_type(file)
            current_section = None

            with open(file, "r") as f:
                for line in f:
                    clean_line = line.split("!")[0].split(";")[0].strip()
                    if not clean_line:
                        continue

                    if filetype in ["prm", "str"]:
                        parts = clean_line.split()

                        # Intercept MASS keyword
                        if parts and parts[0] == "MASS":
                            current_section = "atomtypes"

                        if clean_line in charmm_mapping:
                            current_section = charmm_mapping[clean_line]
                            continue

                        if current_section in self.parameters:
                            if len(parts) >= 2:
                                self._convert_and_store(current_section, parts)

    def _convert_and_store(self, section, parts):
        """
        Converts parameter line to GROMACS format and units, and stores it under
        forward and reverse atom-type tuples.

        :param section: The topology section being processed (e.g., 'bonds', 'angles').
        :type section: str
        :param parts: The split string components of the parsed CHARMM parameter line.
        :type parts: list
        :return: None. Updates the internal self.parameters dictionary.
        """

        formatted_line = None
        keys = []

        try:
            if section == "bonds" and len(parts) >= 4:
                b0_gmx = float(parts[3]) / 10
                kb_gmx = float(parts[2]) * 4.184 * 100 * 2
                formatted_line = (
                    f"{parts[0]:<8} {parts[1]:<8} 1    {b0_gmx:<11.5f} {kb_gmx:.1f}"
                )
                keys = [tuple(parts[0:2]), tuple(reversed(parts[0:2]))]

            elif section == "angles" and len(parts) >= 5:
                k_gmx = float(parts[3]) * 4.184 * 2
                theta_gmx = float(parts[4])

                if len(parts) >= 7:
                    kub_gmx = float(parts[5]) * 4.184 * 100 * 2
                    s0_gmx = float(parts[6]) / 10
                else:
                    kub_gmx = 0.0
                    s0_gmx = 0.0

                formatted_line = f"{parts[0]:<8} {parts[1]:<8} {parts[2]:<8} 5  {theta_gmx:<11.5f} {k_gmx:.1f}   {s0_gmx:.5f}  {kub_gmx:.1f}"
                keys = [tuple(parts[0:3]), tuple(reversed(parts[0:3]))]

            elif section == "dihedrals_proper" and len(parts) >= 7:
                k_gmx = float(parts[4]) * 4.184
                formatted_line = f"{parts[0]:<8} {parts[1]:<8} {parts[2]:<8} {parts[3]:<8} 9    {float(parts[6]):.1f}   {k_gmx:.1f}   {int(parts[5])}"
                keys = [tuple(parts[0:4]), tuple(reversed(parts[0:4]))]

            elif section == "dihedrals_improper" and len(parts) >= 7:
                k_gmx = float(parts[4]) * 4.184 * 2
                formatted_line = f"{parts[0]:<8} {parts[1]:<8} {parts[2]:<8} {parts[3]:<8} 2      {float(parts[6]):.1f}   {k_gmx:.1f}"
                keys = [tuple(parts[0:4]), tuple(reversed(parts[0:4]))]

            elif section == "atomtypes" and len(parts) >= 3 and parts[0] == "MASS":
                atom_type = parts[2]
                keys = [(atom_type,)]
                formatted_line = " ".join(parts)

        except ValueError:
            return

        if formatted_line:
            formatted_line += "\n"
            for k in keys:
                if k not in self.parameters[section]:
                    self.parameters[section][k] = []
                self.parameters[section][k].append(formatted_line)

    def extract_junction_params(self, requested_params: dict):
        """
        Takes a dict of required atom signatures and scrapes them from the database.

        :param requested_params: Dictionary of topology sections mapping to sets of atom type tuples.
        :type requested_params: dict
        :return: None. Populates the internal self.junction_output dictionary.
        """

        # Aliases allow cross-forcefield boundary atoms to find their base parameters
        aliases = {"CC2O1": "CC", "NC2D1": "NH2", "OC2D1": "O", "HCP1": "H"}

        for section, targets in requested_params.items():
            if section not in self.parameters:
                continue

            search_buckets = [section]
            if "dihedral" in section:
                search_buckets = ["dihedrals_proper", "dihedrals_improper"]

            for target_tuple in targets:
                found = False
                fallback_tuple = tuple(aliases.get(atom, atom) for atom in target_tuple)

                for bucket in search_buckets:
                    if found:
                        break

                    for query in [target_tuple, fallback_tuple]:
                        if query in self.parameters[bucket]:
                            for raw_line in self.parameters[bucket][query]:
                                parts = raw_line.split()
                                for i in range(len(target_tuple)):
                                    parts[i] = target_tuple[i]
                                padded_atoms = [
                                    f"{p:<8}" for p in parts[: len(target_tuple)]
                                ]
                                injected_line = (
                                    " ".join(padded_atoms + parts[len(target_tuple) :])
                                    + "\n"
                                )
                                self.junction_output[section].add(injected_line)
                            found = True
                            break

                    if not found and len(target_tuple) == 4:
                        wildcards = [
                            ("X", target_tuple[1], target_tuple[2], "X"),
                            ("X", target_tuple[2], target_tuple[1], "X"),
                        ]
                        for pattern in wildcards:
                            if pattern in self.parameters[bucket]:
                                for raw_line in self.parameters[bucket][pattern]:
                                    parts = raw_line.split()
                                    parts[0] = target_tuple[0]
                                    parts[3] = target_tuple[3]
                                    injected_line = f"{parts[0]:<8} {parts[1]:<8} {parts[2]:<8} {parts[3]:<8} {' '.join(parts[4:])}\n"
                                    self.junction_output[section].add(injected_line)
                                found = True
                                break

                if not found:
                    print(f"WARNING: Missing parameter for {section}: {target_tuple}")

    def write_ff(self, outfile):
        """
        Writes scraped junction parameters to a GROMACS itp file.

        :param outfile: The target filepath to write the topology parameters.
        :type outfile: str
        :return: None. Writes directly to disk.
        """

        # Translator map for strict GROMACS global directives
        gmx_headers = {
            "bonds": "bondtypes",
            "angles": "angletypes",
            "dihedrals_proper": "dihedraltypes",
            "dihedrals_improper": "dihedraltypes",
        }

        grouped_output = {}
        for section, lines in self.junction_output.items():
            if not lines or section not in gmx_headers:
                continue

            header = gmx_headers[section]
            if header not in grouped_output:
                grouped_output[header] = set()
            grouped_output[header].update(lines)

        with open(outfile, "w") as f:
            for header, lines in grouped_output.items():
                f.write(f"[ {header} ]\n")
                for line in sorted(list(lines)):
                    f.write(line)
                f.write("\n")
