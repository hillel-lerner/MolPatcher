import os
import sys
import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from mol_patcher.pdb_io import PdbParser
from mol_patcher.topology_io import TopologyParser
from mol_patcher.mol_record import Mol
from mol_patcher.geometry import MolGraph
from mol_patcher.utilities import get_distance

# --- NEW FUNCTION WRAPPER ---
def compare_matrices(pdb_path, itp_path):
    print(f"\nTesting files: {os.path.basename(pdb_path)} & {os.path.basename(itp_path)}")
    
    pdb_parser = PdbParser()
    _, records, _ = pdb_parser.read_file(pdb_path)
    
    test_mol = Mol("test_mol", records=records)
    test_mol.load_itp(itp_path)

    graph_engine = MolGraph(test_mol)

    matrix_old = graph_engine.matrix.copy()

    graph_engine.gen_topology_matrix()
    matrix_new = graph_engine.matrix.copy()

    if np.array_equal(matrix_old, matrix_new):
        print(" -> PASS")
    else: 
        print(" -> FAIL")
        
        # Find exactly where the arrays disagree
        differences = np.where(matrix_old != matrix_new)
        
        print("    Mismatched Bonds Detected:")
        printed_pairs = set()
        for row, col in zip(differences[0], differences[1]):
            # Use a sorted tuple to avoid printing both directions of the symmetric bond
            pair = tuple(sorted([row, col]))
            if pair not in printed_pairs:
                printed_pairs.add(pair)
                
                atom1 = test_mol.records[pair[0]]
                atom2 = test_mol.records[pair[1]]
                dist = get_distance((atom1.x, atom1.y, atom1.z), (atom2.x, atom2.y, atom2.z))
                
                # Determine which method drew the bond
                if matrix_old[row, col] == 1:
                    status = "Old Method Hallucinated Bond (False Positive)"
                else:
                    status = "Old Method Missed Bond (False Negative)"
                    
                print(f"    [!] {atom1.res_name.strip()}-{atom1.res_seq} {atom1.name.strip()} <-> "
                    f"{atom2.res_name.strip()}-{atom2.res_seq} {atom2.name.strip()} | "
                    f"Dist: {dist:.3f} Å | {status}")
# ----------------------------

if __name__ == "__main__":
    print("Testing MolGraph ITP Connectivity Builder (V2)...")

    # Paths for the original patch
    patch_pdb = os.path.join(project_root, 'pdbs', 'pfp_patch_new.pdb')
    patch_itp = os.path.join(project_root, 'topology_files', 'pfp_patch_new.itp')

    # Paths for the stitched output
    # (Update 'patched_lys_212' to whatever directory your test files are in)
    stitched_pdb = os.path.join(project_root, 'pdbs', 'patched_lys_188.pdb')
    stitched_itp = os.path.join(project_root, 'topology_files', 'patched_lys_188.itp')

    # Run both tests
    compare_matrices(patch_pdb, patch_itp)
    compare_matrices(stitched_pdb, stitched_itp)