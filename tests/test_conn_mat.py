import os
import sys
import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from mol_patcher.pdb_io import PdbParser
from mol_patcher.mol_record import Mol
from mol_patcher.geometry import MolGraph
from mol_patcher.utilities import get_distance

if __name__ == "__main__":
    print("Testing MolGraph Connectivity Builder...\n")

    pdb_dir = os.path.join(project_root, 'pdbs')
    patch_pdb = os.path.join(pdb_dir, "pfp_patch_new.pdb")

    _, records, _, _ = PdbParser.read_file(patch_pdb)
    
    test_mol = Mol("patch", records=records)

    # 2. Build the geometry graph
    graph_engine = MolGraph(test_mol)


    print(f"Total Atoms: {len(test_mol.records)}")
    
    # The matrix is symmetric (A bonded to B means B bonded to A). 
    # Dividing the sum of the matrix by 2 gives the absolute number of unique bonds.
    total_bonds = int(np.sum(graph_engine.matrix) / 2)
    print(f"Total Bonds Found: {total_bonds}")
    
    # If the patch is a single intact molecule, Nmols should be exactly 1.
    print(f"Connected Components: {graph_engine.Nmols}\n")

    # 4. Inspect the chemistry of the bridge atom (C7)
    target_name = "C7"
    
    # Find the integer index of C7 in the PDB records list
    target_idx = next(i for i, r in enumerate(test_mol.records) if r.name.strip() == target_name)
    
    print(f"Checking bonds for {target_name} (List Index {target_idx}):")
    
    # Use np.where to pull all column indices where the matrix has a '1' in C7's row
    bonded_indices = np.where(graph_engine.matrix[target_idx, :] == 1)[0]
    
    target_coords = [test_mol.records[target_idx].x, test_mol.records[target_idx].y, test_mol.records[target_idx].z]

    for b_idx in bonded_indices:
        neighbor = test_mol.records[b_idx]
        neighbor_coords = [neighbor.x, neighbor.y, neighbor.z]
        
        # Calculate the exact distance to verify it fell under the 1.65 / 1.15 cutoffs
        dist = get_distance(target_coords, neighbor_coords)
        print(f"  - Bonded to: {neighbor.name.strip():<4s} | Distance: {dist:.3f} Å")