"""
Handles reading, parsing, formatting, and writing of PDB coordinate files.
"""

import os
from .mol_record import PdbRecord


class PdbParser:
    """
    Handles reading, parsing, and cleaning of PDB coordinate files.
    """

    def fix_chain_id(self, records):
        """
        Identifies missing chain IDs and attempts to patch them using the segment ID.

        :param records: PdbRecord objects parsed from the file.
        :type records: list
        :return: The updated PdbRecord objects.
        :rtype: list
        """

        if not records:
            return records
        needs_fix = any(atom.chain == " " for atom in records)
        if not needs_fix:
            return records

        print("   -> Missing chain IDs detected in PDB. Patching via seg_id...")

        for atom in records:
            chain = atom.chain
            if chain == " ":
                if atom.seg_id:
                    chain = atom.seg_id[-1]
                else:
                    print(
                        f"No chain identified for atom {atom.serial}, setting to blank."
                    )
            atom.chain = chain
        return records

    def fix_res_num(self, records):
        """
        Detects non-sequential residue numbering or insertion codes and dynamically
        renumbers the entire chain sequentially starting from 1.

        :param records: PdbRecord objects to check and renumber.
        :type records: list
        :return: The sequentially renumbered PdbRecord objects.
        :rtype: list
        """

        if not records:
            return records

        needs_fix = False
        expected_counter = 1
        last_seen = (records[0].chain, records[0].res_seq, records[0].ins_code)

        for atom in records:
            current_res = (atom.chain, atom.res_seq, atom.ins_code)
            if current_res[0] != last_seen[0]:  # New chain detected
                expected_counter = 1
                last_seen = current_res
            elif current_res != last_seen:  # New residue detected
                expected_counter += 1
                last_seen = current_res

            # If the actual data doesn't match expected counter, flag it.
            if atom.res_seq != expected_counter or atom.ins_code != " ":
                needs_fix = True
                break

        if not needs_fix:
            return records

        print(
            "   -> Non-sequential numbering or insertion codes detected in PDB. Fixing..."
        )

        self.res_data = []
        self.res_counter = 1
        last_seen = (records[0].chain, records[0].res_seq, records[0].ins_code)

        for atom in records:
            current_res = (atom.chain, atom.res_seq, atom.ins_code)
            if current_res[0] != last_seen[0]:
                self.res_counter = 1
                last_seen = current_res
            elif current_res != last_seen:
                self.res_counter += 1
                last_seen = current_res
            atom.res_seq = self.res_counter
            atom.ins_code = " "
        return records

    def read_file(self, file):
        """
        Reads a PDB file and extracts headers and ATOM/HETATM records.

        :param file: Path to the PDB file.
        :type file: str
        :return: A tuple containing the header lines, PdbRecord objects, and the filename.
        :rtype: tuple of (list, list, str)
        """

        self.file_name = os.path.basename(file)

        self.header_lines = []
        self.atoms = []

        headers = (
            "HEADER",
            "OBSLTE",
            "TITLE",
            "SPLIT",
            "CAVEAT",
            "COMPND",
            "SOURCE",
            "KEYWDS",
            "EXPDTA",
            "NUMMDL",
            "MDLTYP",
            "AUTHOR",
            "REVDAT",
            "SPRSDE",
            "JRNL",
            "REMARK",
            "DBREF",
        )

        with open(file, "r") as f:
            for line in f:
                if line.startswith(headers):
                    self.header_lines.append(line)

                elif line.startswith(("ATOM", "HETATM")):
                    atom = self.parse_line(line, self.file_name)
                    if atom is not None:
                        self.atoms.append(atom)

        return self.header_lines, self.atoms, self.file_name

    def parse_line(self, line, file_name):
        """
        Parses a single PDB line of ATOM or HETATM information into a PdbRecord object.

        :param line: The raw text line from the PDB file.
        :type line: str
        :param file_name: The name of the file being parsed.
        :type file_name: str
        :return: The populated PdbRecord object, or None if parsing fails.
        :rtype: PdbRecord or None
        """

        try:
            return PdbRecord(
                pdb_name=str(file_name),
                record_type=line[0:6].strip(),
                serial=int(line[6:11].strip()),
                name=line[12:16].strip(),
                res_name=line[17:21].strip(),
                chain=line[21],
                res_seq=int(line[22:26].strip()),
                ins_code=line[26],
                x=float(line[30:38]),
                y=float(line[38:46]),
                z=float(line[46:54]),
                seg_id=line[72:76].strip() if len(line) >= 76 else "",
            )

        except ValueError:
            return None


class PdbBuilder:
    """
    Handles formatting and writing PdbRecord objects back into standard PDB file structure.
    """

    def __init__(self, new_pdb_name, atom_list, headers=None):
        """
        Constructs the PdbBuilder.

        :param new_pdb_name: The output file path for the generated PDB.
        :type new_pdb_name: str
        :param atom_list: PdbRecord objects to be written to the file.
        :type atom_list: list
        :param headers: Original PDB header strings to preserve, defaults to None.
        :type headers: list, optional
        """

        self.new_pdb_name = new_pdb_name
        self.atom_list = atom_list
        self.headers = headers if headers else []

    def format_lines(self):
        """
        Formats the internal PdbRecord objects into standard WWPDB text lines.

        :return: A list of formatted PDB text lines ready for writing.
        :rtype: list
        """

        formatted_lines = []

        for header in self.headers:
            formatted_lines.append(header)

        # =========================================================================
        # PDB ATOM RECORD FORMATTING EXPLANATION
        # Reference: https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
        #
        # Col 1-6   : Record Type  (e.g., "ATOM  ")
        # Col 7-11  : Serial No.   (Integer, Right align)
        # Col 12    : Empty        (Handled by the first space in "  ")
        # Col 13-16 : Atom Name    (Centered 4-char field)
        #             * "  {:<3s}" forces name to start at Col 14.
        # Col 17    : AltLoc       (Char)
        # Col 18-20 : ResName      (e.g. "LYS"). NOTE THAT IN CASES OF 4 LETTER RES NAMES (e.g., BGLC), ALTLOC WILL BE REPURPOSED.
        # Col 22    : Chain        (Char)
        # Col 23-26 : ResSeq       (Integer, Right align)
        # Col 27    : InsCode      (Char)
        # Col 31-38 : X Coord      (8.3f)
        # Col 39-46 : Y Coord      (8.3f)
        # Col 47-54 : Z Coord      (8.3f)
        # Col 55-60 : Occupancy    (6.2f)
        # Col 61-66 : TempFactor   (6.2f)
        # Col 73-76 : SegID        (Left align)
        # =========================================================================

        for atom in self.atom_list:
            # PDB format requires occupancy and temp factor even if they are 1.00 and 0.00
            occ = getattr(atom, "occupancy", 1.00)
            temp = getattr(atom, "temp_factor", 0.00)

            line = "{:<6s}{:>5d} {:^4s}{:>4s} {:1s}{:>4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:<4s}\n".format(
                str(atom.record_type),
                int(atom.serial),
                str(atom.name),
                str(atom.res_name)[
                    :4
                ],  # Merged AltLoc and ResName into a single 4-character block {:>4s}. Truncates at 4 so columns never shift.
                str(atom.chain),
                int(atom.res_seq),
                str(atom.ins_code),
                float(atom.x),
                float(atom.y),
                float(atom.z),
                float(occ),
                float(temp),
                str(atom.seg_id),
            )
            formatted_lines.append(line)

        formatted_lines.append("TER\n")
        formatted_lines.append("END\n")

        return formatted_lines

    def write_pdb(self):
        """
        Writes the formatted PDB lines to disk.

        :return: None.
        """

        lines = self.format_lines()
        with open(self.new_pdb_name, "w") as file:
            file.writelines(lines)
        print(f"Patched {self.new_pdb_name}")

