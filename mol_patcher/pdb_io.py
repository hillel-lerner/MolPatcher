import numpy as np
import os
from .mol_record import PdbRecord

class PdbParser:
    """
    Handles reading, parsing, and cleaning of PDB coordinate files.
    """

    def fix_chain_id(self, records):
        """
        Identifies missing chain IDs and attempts to patch them using the segment ID.
        
        :param records: (list) PdbRecord objects parsed from the file.
        :return: (list) The updated PdbRecord objects.
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
                    print(f"No chain identified for atom {atom.serial}, setting to blank.")
            atom.chain = chain
        return records

    def fix_res_num(self, records):
        """
        Detects non-sequential residue numbering or insertion codes and dynamically 
        renumbers the entire chain sequentially starting from 1.
        
        :param records: (list) PdbRecord objects.
        :return: (list) The renumbered PdbRecord objects.
        """

        if not records:
            return records

        needs_fix = False
        expected_counter = 1
        last_seen = (records[0].chain, records[0].res_seq, records[0].ins_code)
        
        for atom in records:
            current_res = (atom.chain, atom.res_seq, atom.ins_code)
            if current_res[0] != last_seen[0]:    # New chain detected
                expected_counter = 1
                last_seen = current_res
            elif current_res != last_seen:        # New residue detected
                expected_counter += 1
                last_seen = current_res
                
            # If the actual data doesn't match expected counter, flag it.
            if atom.res_seq != expected_counter or atom.ins_code != " ":
                needs_fix = True
                break
                
        if not needs_fix:
            return records
            
        print("   -> Non-sequential numbering or insertion codes detected in PDB. Fixing...")
        
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
        
        :param file: (str) Path to the PDB file.
        :return: (tuple) A tuple containing (header_lines_list, PdbRecord_list, file_name_string).
        """

        self.file_name = os.path.basename(file)

        self.header_lines = []
        self.atoms = []

        headers = ('HEADER', 'OBSLTE', 'TITLE', 'SPLIT', 'CAVEAT', 'COMPND', 'SOURCE', 
                'KEYWDS', 'EXPDTA', 'NUMMDL', 'MDLTYP', 'AUTHOR', 'REVDAT', 
                'SPRSDE', 'JRNL', 'REMARK', 'DBREF')

        with open(file, 'r') as f:
            for line in f:
                if line.startswith(headers):
                    self.header_lines.append(line)

                elif line.startswith(('ATOM', 'HETATM')):
                    atom = self.parse_line(line, self.file_name) 
                    if atom is not None:
                        self.atoms.append(atom)

        return self.header_lines, self.atoms, self.file_name


    def parse_line(self, line, file_name):

        """parses pdb line of ATOM or HETATM information into an PdbRecord object"""
        try:
            return PdbRecord(
                pdb_name=str(file_name),
                record_type=line[0:6].strip(),
                serial=int(line[6:11].strip()),
                name=line[12:16].strip(),
                res_name=line[16:20].strip(),
                chain=line[21],
                res_seq=int(line[22:26].strip()),
                ins_code=line[26],
                x=float(line[30:38]),
                y=float(line[38:46]),
                z=float(line[46:54]),
                seg_id=line[72:76].strip() if len(line) >= 76 else ""
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
        
        :param new_pdb_name: (str) The output file path.
        :param atom_list: (list) PdbRecord objects to be written.
        :param headers: (list, optional) Original PDB header strings to preserve.
        """
        
        self.new_pdb_name = new_pdb_name
        self.atom_list = atom_list
        self.headers = headers if headers else []

    def format_lines(self):
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
            occ = getattr(atom, 'occupancy', 1.00)
            temp = getattr(atom, 'temp_factor', 0.00)
            
            line = "{:<6s}{:>5d} {:^4s}{:>4s} {:1s}{:>4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:<4s}\n".format(
                str(atom.record_type),
                int(atom.serial),   
                str(atom.name),
                str(atom.res_name)[:4], # Merged AltLoc and ResName into a single 4-character block {:>4s}. Truncates at 4 so columns never shift. 
                str(atom.chain),
                int(atom.res_seq),  
                str(atom.ins_code), 
                float(atom.x), float(atom.y), float(atom.z),
                float(occ),
                float(temp),
                str(atom.seg_id)
            )
            formatted_lines.append(line)

        formatted_lines.append("TER\n")
        formatted_lines.append("END\n")

        return formatted_lines

    def write_pdb(self):
        lines = self.format_lines()
        with open(self.new_pdb_name, 'w') as file:
            file.writelines(lines)
        print(f"Patched {self.new_pdb_name}")