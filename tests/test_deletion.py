import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from mol_patcher.pdb_io import PdbParser
from mol_patcher.mol_record import Mol
from mol_patcher.geometry import MolGraph

if __name__ == "__main__":
    print("Testing Graph-Based Atom Deletion...\n")

    # Load the protein base
    pdb_dir = os.path.join(project_root, 'pdbs')
    base_pdb = os.path.join(pdb_dir, "step3_input.pdb")

    _, records, _, _ = PdbParser.read_file(base_pdb)
    base_mol = Mol("protein", records=records)

    # Build the geometry graph
    print("Building full protein graph (this may take a few seconds)...")
    base_graph = MolGraph(base_mol)

    # Find the NZ anchor on the target residue
    target_res = 136
    anchor_name = "NZ"
    
    # Get the exact list index of the anchor atom
    anchor_idx = next(i for i, r in enumerate(base_mol.records) 
                    if r.res_seq == target_res and r.name.strip() == anchor_name)
    
    print(f"\nAnalyzing neighbors for {anchor_name} on Residue {target_res}:")

    # 4. Use the graph to instantly fetch all bonded atoms
    neighbors = base_graph.nx_graph.neighbors(anchor_idx)
    
    # 5. Filter the neighbors to find only the Hydrogens
    hydrogens_to_delete = []
    for n_idx in neighbors:
        neighbor_record = base_mol.records[n_idx]
        print(f"  - Found bonded atom: {neighbor_record.name.strip()}")
        
        if neighbor_record.name.strip().startswith("H"):
            hydrogens_to_delete.append(neighbor_record)
            
    print(f"\nGraph identified {len(hydrogens_to_delete)} hydrogens for deletion.")