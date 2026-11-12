"""
Handles parsing, cleaning, and writing of GROMACS .itp topology files.
"""

import os

from .mol_record import ItpAngle, ItpAtom, ItpBond, ItpCmap, ItpDih, ItpPair
from .topology_tools import reindex_topology


class TopologyParser:
    """
    Handles parsing and cleanup of GROMACS .itp topology files.
    """

    @staticmethod
    def clean_itp_file(itp_path, temp):
        """
        Scans the ITP file for non-integer residue designations (e.g., '52A').
        If found, generates a sanitized temporary copy with sequential integers to prevent parser crashes.

        :param itp_path: Path to the input .itp file.
        :type itp_path: str
        :param temp: Path to the output directory for the temporary sanitized file.
        :type temp: str
        :return: Path to the clean .itp file (returns original path if no fix was needed).
        :rtype: str
        """

        needs_fix = False
        with open(itp_path, "r") as f:
            in_atoms = False
            for line in f:
                stripped = line.strip()
                if stripped.startswith("[ atoms ]"):
                    in_atoms = True
                    continue
                elif stripped.startswith("[") and in_atoms:
                    break

                if in_atoms and stripped and not stripped.startswith(";"):
                    parts = stripped.split()
                    if len(parts) > 2 and not parts[2].isdigit():
                        needs_fix = True
                        break

        # If all residues are integers, proceed with the original file
        if not needs_fix:
            return itp_path

        print(
            f"Warning: Non-integer residue detected in {os.path.basename(itp_path)}. Cleaning topology before loading..."
        )
        clean_itp_path = os.path.join(temp, "temp_sanitized.itp")

        with open(itp_path, "r") as infile, open(clean_itp_path, "w") as outfile:
            in_atoms = False
            res_counter = 1
            last_seen_raw = None

            for line in infile:
                if line.strip().startswith("[ atoms ]"):
                    in_atoms = True
                    outfile.write(line)
                    continue
                elif line.strip().startswith("[") and in_atoms:
                    in_atoms = False

                if in_atoms and line.strip() and not line.strip().startswith(";"):
                    parts = line.split()
                    raw_res = parts[2]

                    base_int_str = "".join(filter(str.isdigit, raw_res))
                    base_int = int(base_int_str) if base_int_str else 0

                    if last_seen_raw is None:
                        last_seen_raw = raw_res
                        res_counter = 1
                    elif base_int == 1 and raw_res != last_seen_raw:
                        res_counter = 1
                        last_seen_raw = raw_res
                    elif raw_res != last_seen_raw:
                        res_counter += 1
                        last_seen_raw = raw_res

                    # Reconstruct the line
                    comment_append = ""
                    if str(res_counter) != raw_res:
                        comment_append = f" ; old_res: {raw_res}"

                    formatted_line = f"{parts[0]:>6} {parts[1]:>10} {res_counter!s:>6} {parts[3]:>6} {parts[4]:>6} {parts[5]:>6} {parts[6]:>12} {parts[7]:>12}{comment_append}"
                    outfile.write(formatted_line + "\n")
                else:
                    outfile.write(line)

        return clean_itp_path

    @staticmethod
    def read_file(filepath):
        """
        Reads a clean GROMACS .itp file and extracts topological parameters.

        :param filepath: Path to the .itp file.
        :type filepath: str
        :return: A tuple containing the moleculetype string and lists of topology dataclasses.
        :rtype: tuple of (str, list, list, list, list, list, list)
        """

        atoms, bonds, pairs, angles, dihs, cmap = [], [], [], [], [], []

        section_map = {
            "bonds": (bonds, ItpBond, 3),
            "pairs": (pairs, ItpPair, 3),
            "angles": (angles, ItpAngle, 4),
            "dihedrals": (dihs, ItpDih, 5),
            "cmap": (cmap, ItpCmap, 6),
        }

        current_section = None
        moltype_lines = []

        with open(filepath, "r") as f:
            for line in f:
                stripped = line.strip()

                if stripped.startswith("["):
                    current_section = stripped.strip("[] ").lower()
                    if current_section == "moleculetype":
                        moltype_lines.append(line)
                    continue

                if current_section == "moleculetype":
                    moltype_lines.append(line)
                    continue

                # Ignore blank lines and comments
                if not stripped or stripped.startswith(";"):
                    continue

                parts = line.split()

                if (
                    current_section == "atoms"
                    and len(parts) >= 8
                    and parts[0].isdigit()
                ):
                    comment = ""
                    if ";" in line:
                        raw_comment = line[line.find(";") :].strip()
                        if "old_res" in raw_comment:
                            comment = (
                                "; old_res: "
                                + raw_comment.split("old_res:")[-1].strip()
                            )

                    atoms.append(
                        ItpAtom(
                            int(parts[0]),
                            parts[1],
                            int(parts[2]),
                            parts[3],
                            parts[4],
                            int(parts[5]),
                            float(parts[6]),
                            float(parts[7]),
                            comment,
                        )
                    )

                elif current_section in section_map:
                    target_list, data_class, required_len = section_map[current_section]
                    if len(parts) >= required_len:
                        target_list.append(
                            data_class(*[int(p) for p in parts[:required_len]])
                        )

            moltype_section = "".join(moltype_lines)
            return moltype_section, atoms, bonds, pairs, angles, dihs, cmap


class TopologyBuilder:
    """
    Handles formatting and writing topology objects back into a GROMACS-compatible .itp file.
    """

    def __init__(self, mol, filename, mol_name):
        """
        Constructs the TopologyBuilder.

        :param mol: The molecule object containing the topology to write.
        :type mol: Mol
        :param filename: The output file path for the .itp file.
        :type filename: str
        :param mol_name: The molecule name to be written in the [ moleculetype ] header.
        :type mol_name: str
        """

        self.mol = mol
        self.filename = filename
        self.mol_name = mol_name

    def write_itp(self, posre_fc=1000):
        """
        Synchronizes internal atom indices and writes the topology out to disk.
        :param posre_fc: The force constant for the position restraints section of the itp in kJ/mol/nm2

        :return: None.
        """

        reindex_topology(self.mol)

        with open(self.filename, "w") as f:
            # Notes
            f.writelines(f"; {note}\n" for note in self.mol.notes)
            if self.mol.notes:
                f.write("\n")

            # Molecule Header
            if self.mol.moltype_section:
                f.write(self.mol.moltype_section)
                if not self.mol.moltype_section.endswith("\n\n"):
                    f.write("\n")
            else:
                f.write("[ moleculetype ]\n; Name            nrexcl\n")
                f.write(f"{self.mol.name:<15} 3\n\n")

            # Atoms Section
            if self.mol.atoms:
                f.write("[ atoms ]\n; nr type resnr residue atom cgnr charge mass\n")
                for atom in self.mol.atoms:
                    # Check if the atom object has a comment, default to empty string if not
                    comment = getattr(atom, "comment", "")
                    f.write(
                        f"{atom.number:>6d} {atom.type.strip():>10s} {atom.res_n:>6d} {atom.res:>6s} "
                        f"{atom.atom:>6s} {atom.cgnr:>6d} {atom.charge:>10.4f} {atom.mass:>10.4f} {comment}\n"
                    )
                f.write("\n")

            # Bonds Section
            if self.mol.bonds:
                f.write("[ bonds ]\n;  ai    aj  funct\n")

                f.writelines(
                    f"{b.a1:>7d} {b.a2:>7d} {b.type:>7d}\n" for b in self.mol.bonds
                )
                f.write("\n")

            # Pairs Section
            if self.mol.pairs:
                f.write("[ pairs ]\n;  ai    aj  funct\n")
                f.writelines(
                    f"{p.a1:>7d} {p.a2:>7d} {p.type:>7d}\n" for p in self.mol.pairs
                )
                f.write("\n")

            # Angles Section
            if self.mol.angles:
                f.write("[ angles ]\n;  ai    aj    ak  funct\n")
                f.writelines(
                    f"{a.a1:>7d} {a.a2:>7d} {a.a3:>7d} {a.type:>7d}\n"
                    for a in self.mol.angles
                )
                f.write("\n")

            # Dihedrals Section
            if self.mol.dihs:
                f.write("[ dihedrals ]\n;  ai    aj    ak    al  funct\n")
                f.writelines(
                    f"{d.a1:>7d} {d.a2:>7d} {d.a3:>7d} {d.a4:>7d} {d.type:>7d}\n"
                    for d in self.mol.dihs
                )
                f.write("\n")

            # CMAP section
            if self.mol.cmap:
                f.write("[ cmap ]\n;  ai    aj    ak    al  am funct\n")
                f.writelines(
                    f"{cmp.a1:>7d} {cmp.a2:>7d} {cmp.a3:>7d} {cmp.a4:>7d} {cmp.a5:>7d} {cmp.type:>7d}\n"
                    for cmp in self.mol.cmap
                )
                f.write("\n")

            f.write(
                "#ifdef POSRES\n[ position_restraints ]\n;  ai  funct   fcx    fcy    fcz\n"
            )
            f.writelines(
                f"{a.number:>6d}    1 {posre_fc:>7} {posre_fc:>7} {posre_fc:>7}\n"
                for a in self.mol.atoms
                if a.mass > 1.5
            )
            f.write("#endif\n")


def top_lines(ff_name, includes, mol_name, system_name):
    """
    Builds a list of lines for the topol.top file needed for gromacs runs with the patched itp files.

    :param ff_name: Path to the forcefield file, if it exists.
    :type ff_name: str
    :param includes: List of paths to the forcefield and topology includes relative to the .top file.
    :type includes: list[str]
    :param mol_name: The name of the molecule listed in the molecule section of topol.top.
    :type mol_name: str
    :param system_name: The name of the system listed in the system section of topol.top.
    :type system_name: str
    :return top: A list of lines that will be written to `topol.top`
    :rtype: list
    """

    lines = []

    lines.append(";;")
    lines.append(";; Generated by MolPatcher for the patched molecule only.")
    lines.append(";; Solvent, ions and any other molecules must be added by hand:")
    lines.append(";;  add their #include lines below the patched molecule's include;")
    lines.append(
        ";;  gmx solvate -p and gmx genion -p fill in their [ molecules ] counts."
    )
    lines.append(";;")
    lines.append(";; The main GROMACS topology file")
    lines.append(";;")
    lines.append("")
    lines.append("; Include forcefield parameters")
    if ff_name:
        if ff_name.endswith(".itp"):
            ff_include = ff_name
        elif ff_name.endswith(".ff"):
            ff_include = f"{ff_name}/forcefield.itp"
        else:
            ff_include = f"{ff_name}.ff/forcefield.itp"
        lines.append(f'#include "{ff_include}"')
    else:
        lines.append(
            '; #include "<forcefield>.ff/forcefield.itp"   <- set your force field'
        )
    lines.extend(f'#include "{i}"' for i in includes)
    lines.append("")
    lines.append("[ system ]")
    lines.append(f"{system_name}")
    lines.append("")
    lines.append("[ molecules ]")
    lines.append(f"{mol_name:<15} 1")

    return lines


def write_top(path, ff_name, includes, mol_name, system_name):
    """
    Writes the lines build in top_lines to a file

    :param path: Path to topol.top
    :type path: str
    :param ff_name: Path to the forcefield file, if it exists.
    :type ff_name: str
    :param includes: List of paths to the forcefield and topology includes relative to the .top file.
    :type includes: list[str]
    :param mol_name: The name of the molecule listed in the molecule section of topol.top.
    :type mol_name: str
    :param system_name: The name of the system listed in the system section of topol.top.
    :type system_name: str
    """

    lines = top_lines(ff_name, includes, mol_name, system_name)

    with open(path, "w") as file:
        file.write("\n".join(lines) + "\n")
