import sys
import os

current_script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_script_dir)

if project_root not in sys.path:
    sys.path.append(project_root)
itp_dir = os.path.join(project_root, 'itps')

from mol_patcher.itp_io import *
from mol_patcher.mol_record import Mol 

def main():
    print("--- Running Test: ITP I/O Identity ---")
    
    infile = os.path.join(itp_dir, "6oge_prob.itp") # Ensure this file exists
    outfile = os.path.join(itp_dir, "6oge_test_io.itp")

    print(f"Reading {infile}...")
    
   
    atoms, bonds, pairs, angles, dihs = ItpParser.read_file(infile)

    mol_obj = Mol("6oge_prob")
    mol_obj.atoms = atoms
    mol_obj.bonds = bonds
    mol_obj.pairs = pairs
    mol_obj.angles = angles
    mol_obj.dihs = dihs

    print(f"  - Parsed {len(mol_obj.atoms)} atoms and {len(mol_obj.bonds)} bonds")

    print(f"Writing {outfile}...")
   
   
    builder = BuildItp(mol_obj, outfile)
    

    builder.write_itp() 
    print("  - Identity write completed.")

if __name__ == "__main__":
    main()