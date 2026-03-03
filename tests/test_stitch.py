import os
import sys
import numpy as np
from scipy.spatial.transform import Rotation
from dataclasses import replace

current_test_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_test_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

pdb_dir = os.path.join(project_root, 'pdbs')

from mol_patcher.pdb_io import PdbParser, BuildPdb
from mol_patcher.align_geom import PatchAligner

def main():
    print("--- Running Test: Stitching Gauche Butane ---")
    
    patch_infile = os.path.join(pdb_dir, "butane_01.pdb")  # Tail (to be rotated)
    base_infile  = os.path.join(pdb_dir, "butane_02.pdb")  # Head (static reference)
    outfile      = os.path.join(pdb_dir, "butane_gauche.pdb")
    
    # FIX: Catch all 4 return values from PdbParser.read_file
    headers, patch_atoms, ters, _ = PdbParser.read_file(patch_infile)
    _, base_atoms, _, _ = PdbParser.read_file(base_infile)

    # 1. Define Anchors for Alignment
    # We use C2 and C3 as the common axis for both molecules
    p_anchors = [
        next(a for a in patch_atoms if a.name == "C2"),
        next(a for a in patch_atoms if a.name == "C3"),
        next(a for a in patch_atoms if a.name == "C1") # Third point to lock plane
    ]
    b_anchors = [
        next(a for a in base_atoms if a.name == "C2"),
        next(a for a in base_atoms if a.name == "C3"),
        next(a for a in base_atoms if a.name == "C1")
    ]

    # 2. Align patch to base first
    aligner = PatchAligner(patch_atoms, p_anchors, b_anchors)
    patch_atoms = aligner.implement_align()

    # 3. Define Pivot Axis (C2-C3) from the base
    c2 = next(a for a in base_atoms if a.name == "C2")
    c3 = next(a for a in base_atoms if a.name == "C3")
    
    pivot_point = np.array([c2.x, c2.y, c2.z])
    axis_vector = np.array([c3.x - c2.x, c3.y - c2.y, c3.z - c2.z])
    axis_vector = axis_vector / np.linalg.norm(axis_vector) 
    
    # 4. Create 60-degree rotation matrix
    twist = Rotation.from_rotvec(axis_vector * np.radians(60.0))

    # Apply rotation to patch coordinates
    patch_coords = np.array([[a.x, a.y, a.z] for a in patch_atoms])
    patch_coords = patch_coords - pivot_point     
    patch_coords = twist.apply(patch_coords)      
    patch_coords = patch_coords + pivot_point     

    # Update Patch Objects using 'replace'
    updated_patch = []
    for i, atom in enumerate(patch_atoms):
        updated_patch.append(replace(atom, x=patch_coords[i][0], y=patch_coords[i][1], z=patch_coords[i][2]))

    # 5. STITCH (Combine Halves)
    final_atoms = []
    head_names = ["C1", "H1", "H2", "H3", "C2", "H4", "H5"]
    tail_names = ["C3", "H6", "H7", "C4", "H8", "H9", "H10"]

    for atom in base_atoms:
        if atom.name in head_names:
            final_atoms.append(atom)
            
    for atom in updated_patch:
        if atom.name in tail_names:
            final_atoms.append(atom)

    # 6. Build and Save
    BuildPdb(outfile, final_atoms, headers, ters).write_pdb()
    print(f"Done. Gauche conformer saved to {outfile}")

if __name__ == "__main__":
    main()