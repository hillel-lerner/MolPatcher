import os


class ForceField:
    def __init__(self, ff_files: list) -> None:
        """
        Accepts a list of file paths for multiple CHARMM parameter files to be loaded into a searchable database.
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

    def find_ff_type(self, ff_file):
        base_name = os.path.basename(ff_file)
        filetype = None

        if ".prm" in base_name:
            filetype = "prm"
        elif ".str" in base_name:
            filetype = "str"
        elif ".itp" in base_name:
            filetype = "itp"

        return filetype

    def read_database(self):
        """
        Loops through all provided files and builds the searchable parameter dictionaries.
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
                        if clean_line in charmm_mapping:
                            current_section = charmm_mapping[clean_line]
                            continue

                        if current_section in self.parameters:
                            parts = clean_line.split()

                            if len(parts) >= 2:
                                self._convert_and_store(current_section, parts)

    def _convert_and_store(self, section, parts):
        """
        Converts parameter line to GROMACS format + units and stores it under forward and reverse atom-type tuples
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
                # CHARMM Format: atom1 atom2 atom3 Ktheta Theta0
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
                # CHARMM Format: atom1 atom2 atom3 atom4 Kchi n delta(phase)
                k_gmx = float(parts[4]) * 4.184
                formatted_line = f"{parts[0]:<8} {parts[1]:<8} {parts[2]:<8} {parts[3]:<8} 9    {float(parts[6]):.1f}   {k_gmx:.1f}   {int(parts[5])}"
                keys = [tuple(parts[0:4]), tuple(reversed(parts[0:4]))]

            elif section == "dihedrals_improper" and len(parts) >= 7:
                # CHARMM Format: atom1 atom2 atom3 atom4 Kchi n delta(phase)
                k_gmx = float(parts[4]) * 4.184 * 2
                formatted_line = f"{parts[0]:<8} {parts[1]:<8} {parts[2]:<8} {parts[3]:<8} 2      {float(parts[6]):.1f}   {k_gmx:.1f}"
                keys = [tuple(parts[0:4]), tuple(reversed(parts[0:4]))]

            elif section == "atomtypes" and len(parts) >= 3 and parts[0] == "MASS":
                # CHARMM Format: MASS code atom_type mass
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
        :param requested_params: e.g. {"angles": [("CT2", "CT2", "NG311")], "dihedrals_proper": [...]}
        """
        for section, targets in requested_params.items():
            if section not in self.parameters:
                continue
            for target_tuple in targets:
                if target_tuple in self.parameters[section]:
                    for line in self.parameters[section][target_tuple]:
                        self.junction_output[section].add(line)
                    continue

                found = False

                if (
                    section in ["dihedrals_proper", "dihedrals_improper"]
                    and len(target_tuple) == 4
                ):
                    a1, a2, a3, a4 = target_tuple
                    wildcard_patterns = [("X", a2, a3, "X"), ("X", a3, a2, "X")]

                    for pattern in wildcard_patterns:
                        if pattern in self.parameters[section]:
                            for raw_line in self.parameters[section][pattern]:
                                parts = raw_line.split()
                                parts[0] = a1
                                parts[3] = a4

                                injected_line = f"{parts[0]:<8} {parts[1]:<8} {parts[2]:<8} {parts[3]:<8} {' '.join(parts[4:])}\n"
                                self.junction_output[section].add(injected_line)
                            found = True
                            break

                if not found:
                    print(f"WARNING: Missing parameter for {section}: {target_tuple}")

    def write_ff(self, outfile):
        """
        Writes scraped junction parameters to a GROMACS itp file.
        """

        with open(outfile, "w") as f:
            for section, lines in self.junction_output.items():
                if not lines:
                    continue

                f.write(f"[ {section} ]\n")

                for line in sorted(list(lines)):
                    f.write(line)
                f.write("\n")


# --- Example Execution ---
# ff = ForceField(["par_all36_prot.prm"])
# ff.read_database()
# requests = {
#     "angles": [("CT2", "CT2", "NG311")],
#     "dihedrals_proper": [("CT2", "NG311", "CG301", "CG321")]
# }
# ff.extract_junction_params(requests)
# ff.write_ff("junction_forcefield.itp")
