import os
import re


class ForceField:
    def __init__(self, ff_file) -> None:
        self.ff_file = ff_file
        self.ff_type, self.filetype = self.find_ff_type(self.ff_file)

        self.parameters = {
            "atomtypes": [],
            "bonds": [],
            "angles": [],
            "dihedrals_proper": [],
            "dihedrals_improper": [],
        }

    def find_ff_type(self, ff_file):
        base_name = os.path.basename(ff_file)
        filetype = None
        ff_type = None

        if ".prm" in base_name:
            filetype = "prm"
        elif ".str" in base_name:
            filetype = "str"
        elif ".itp" in base_name:
            filetype = "itp"

        if base_name == "par_all36_carb.prm":
            ff_type = "charmm36_carb"
        elif base_name == "par_all36_prot.prm":
            ff_type = "charmm36_prot"
        elif base_name == "toppar_all36_carb_glycopeptide.str":
            ff_type = "charmm_glycopeptide"

        return ff_type, filetype

    def read_ff(self):

        current_section = None
        itp_header_pattern = re.compile(r"\[\s*(\w+)\s*\]")

        charmm_mapping = {
            "ATOMS": "atomtypes",
            "BONDS": "bonds",
            "ANGLES": "angles",
            "DIHEDRALS": "dihedrals_proper",
            "IMPROPER": "dihedrals_improper",
        }

        with open(self.ff_file, "r") as f:
            for line in f:
                clean_line = line.split("!")[0].split(";")[0].strip()
                if not clean_line:
                    continue

                if self.filetype == "itp":
                    match = itp_header_pattern.match(clean_line)
                    if match:
                        current_section = match.group(1)
                        continue
                elif self.filetype in ["prm", "str"]:
                    if clean_line in charmm_mapping:
                        current_section = charmm_mapping[clean_line]
                        continue
                if current_section in self.parameters:
                    self.parameters[current_section].append(clean_line)

    def inject_params(self, new_params: dict):
        existing_atoms = set()
        existing_bonds = set()
        existing_angles = set()
        existing_dihedrals = set()
        existing_improper = set()

        for line in self.parameters["dihedrals_proper"]:
            parts = line.split()
            atoms = tuple(parts[0:4])
            existing_dihedrals.add(atoms)

    def write_ff(self, outfile):
        pass
