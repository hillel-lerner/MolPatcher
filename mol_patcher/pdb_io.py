import numpy as np
import os
from .mol_record import PdbRecord

class PdbParser:

    @staticmethod
    def read_file(file):

        file_name = os.path.basename(file)

        header_lines = []
        atoms = []
        ter_lines = [] # necessary if ATOM record exists

        headers = ('HEADER', 'OBSLTE', 'TITLE', 'SPLIT', 'CAVEAT', 'COMPND', 'SOURCE', 
                'KEYWDS', 'EXPDTA', 'NUMMDL', 'MDLTYP', 'AUTHOR', 'REVDAT', 
                'SPRSDE', 'JRNL', 'REMARK', 'DBREF')

        with open(file, 'r') as f:
            for line in f:
                if line.startswith(headers):
                    header_lines.append(line)

                elif line.startswith(('ATOM', 'HETATM')):
                    atom = PdbParser.parse_line(line, file_name) 
                    atoms.append(atom)

                elif line.startswith('TER'):
                    ter_lines.append(line)

        return header_lines, atoms, ter_lines, file_name


    @staticmethod
    def parse_line(line, file_name):

        """parses pdb line of ATOM or HETATM information into an PdbRecord object"""
        try:
            return PdbRecord(
                pdb_name=str(file_name),
                record_type=line[0:6].strip(),
                serial=int(line[6:11].strip()),
                name=line[12:16].strip(),
                res_name=line[17:20].strip(),
                chain=line[21],
                res_seq=int(line[22:26].strip()),
                x=float(line[30:38]),
                y=float(line[38:46]),
                z=float(line[46:54]),
                seg_id=line[72:76].strip() if len(line) >= 76 else ""
                )

        except ValueError:
            return None
        
class BuildPdb:

    def __init__(self, new_pdb_name, atom_list, headers=None, ter_line=None):
        self.new_pdb_name = new_pdb_name
        self.atom_list = atom_list
        self.headers = headers if headers else []
        self.ter_line = ter_line if ter_line else []


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
        # Col 18-20 : ResName      (e.g. "LYS")
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
            
            line = "{:<6s}{:>5d} {:^4s}{:1s}{:>3s} {:1s}{:>4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:<4s}\n".format(
                str(atom.record_type),
                int(atom.serial),   
                str(atom.name),
                " ", # AltLoc
                str(atom.res_name),
                str(atom.chain),
                int(atom.res_seq),  
                " ", # InsCode
                float(atom.x), float(atom.y), float(atom.z),
                float(occ),
                float(temp),
                str(atom.seg_id)
            )
            formatted_lines.append(line)

        for ter in self.ter_line:
            formatted_lines.append(ter)
            
        formatted_lines.append("END\n")

        return formatted_lines

    def write_pdb(self):
        lines = self.format_lines()
        with open(self.new_pdb_name, 'w') as file:
            file.writelines(lines)
        print(f"Patched {self.new_pdb_name}")
