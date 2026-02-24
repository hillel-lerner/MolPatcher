import sys
import os

current_script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_script_dir)

if project_root not in sys.path:
    sys.path.append(project_root)
pdb_dir = os.path.join(project_root, 'pdbs')

from mol_patcher.pdb_io import PdbParser, BuildPdb

def main():
    print("--- Running Test 01: PDB I/O Identity ---")
    
    infile = os.path.join(pdb_dir, "6oge_prob_allatom.pdb") # Ensure this file exists
    outfile = os.path.join(pdb_dir, "6oge_test_io.pdb")

    if not os.path.exists(infile):
        print(f"ERROR: Could not find input file at: {infile}")
        return

    print(f"Reading {infile}...")
    # FIX: Catching all 4 return values (headers, atoms, ters, file_name)
    headers, atoms, ters, _ = PdbParser.read_file(infile)
    print(f"  - Found {len(headers)} header lines")
    print(f"  - Found {len(atoms)} atoms")

    print(f"Writing {outfile}...")
    # BuildPdb expects (filename, atoms, headers, ter_line)
    builder = BuildPdb(outfile, atoms, headers, ters)
    builder.write_pdb()
    print("  - Identity write completed.")

if __name__ == "__main__":
    main()