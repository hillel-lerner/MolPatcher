import os
import sys

current_test_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_test_dir)

if project_root not in sys.path:
    sys.path.append(project_root)

pdb_dir = os.path.join(project_root, 'pdbs')

from mol_patcher.pdb_io import PdbParser, BuildPdb
from mol_patcher.align_geom import PatchAligner 

def main():
    print("--- Running Test: Alignment ---") # If you don't see this, the script isn't even starting
    
    file_orig = os.path.join(pdb_dir, "butane_01.pdb")
    file_target = os.path.join(pdb_dir, "butane_02.pdb")
    file_aligned = os.path.join(pdb_dir, "butane_new.pdb")

    # 1. Read files - Catching all 4 return values
    try:
        headers, atoms_patch, ters, _ = PdbParser.read_file(file_orig) #
        _, atoms_base, _, _ = PdbParser.read_file(file_target) #
    except Exception as e:
        print(f"Error reading PDB files: {e}")
        return

    # 2. Define Anchors 
    # Ensure these names (C1, C2, C3) exist in your butane files
    try:
        mobile_anchors = [
            next(a for a in atoms_patch if a.name == "C1"),
            next(a for a in atoms_patch if a.name == "C2"),
            next(a for a in atoms_patch if a.name == "C3")
        ]
        target_anchors = [
            next(a for a in atoms_base if a.name == "C1"),
            next(a for a in atoms_base if a.name == "C2"),
            next(a for a in atoms_base if a.name == "C3")
        ]
    except StopIteration:
        print("Error: Could not find anchor atoms (C1, C2, C3) in butane files.")
        return

    print(f"Aligning {len(atoms_patch)} atoms to {len(atoms_base)} atoms...")
    
    # 3. Align (Using all 3 required arguments)
    aligner = PatchAligner(atoms_patch, mobile_anchors, target_anchors)
    atoms_aligned = aligner.implement_align()

    # 4. Write result
    print(f"Writing result to {os.path.basename(file_aligned)}...")
    BuildPdb(file_aligned, atoms_aligned, headers, ters).write_pdb()

    print("\n[Verification Success]")

if __name__ == "__main__":
    main()