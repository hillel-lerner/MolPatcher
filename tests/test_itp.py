import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from mol_patcher.mol_record import Mol
from mol_patcher.pdb_io import PdbParser
from mol_patcher import mol_stitcher

def test_full_system_integration():
    test_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(test_dir)
    
    protein_pdb = os.path.join(project_root, 'pdbs', '6oge_prob_allatom.pdb')
    protein_itp = os.path.join(project_root, 'itps', '6oge_prob.itp')
    
    print("\n" + "="*60)
    print("STARTING TOPOLOGY SYNCHRONIZATION TEST")
    print("="*60)

    # Load Protein
    headers, target_atoms, ters, t_name = PdbParser.read_file(protein_pdb)
    base_mol = Mol(name=t_name, records=target_atoms)
    base_mol.load_itp(protein_itp)
    
    print(f"Base Protein: {len(base_mol.records)} PDB atoms | {len(base_mol.atoms)} ITP atoms")

    # Load Patch
    pfp_atoms = mol_stitcher.get_pfp_pdb()
    pfp_mol = Mol(name="PFP", records=pfp_atoms)
    mol_stitcher.get_pfp_itp(pfp_mol)

    # Identify Anchors
    target_anchors = [
        base_mol.get_atom(188, " ", "LYS", "CE"),
        base_mol.get_atom(188, " ", "LYS", "CD"),
        base_mol.get_atom(188, " ", "LYS", "NZ")
    ]
    
    # Execute Surgery
    # Note: We pass pfp_mol to ensure its internal bonds are available for reindexing
    stitched_mol = mol_stitcher.stitch_molecules(
        base_mol, pfp_atoms, target_anchors[0], target_anchors, pfp_mol
    )

    # Validate 1: Synchronized lengths
    pdb_len = len(stitched_mol.records)
    itp_len = len(stitched_mol.atoms)
    print(f"Final System: {pdb_len} PDB atoms | {itp_len} ITP atoms")
    assert pdb_len == itp_len, f"Sync Failure! PDB({pdb_len}) != ITP({itp_len})"

    # Validate 2: Correct Re-indexing (1 to N)
    for i, at in enumerate(stitched_mol.atoms):
        assert at.number == i + 1, f"ITP Indexing break at atom {i+1}"
    
    for i, rec in enumerate(stitched_mol.records):
        assert rec.serial == i + 1, f"PDB Serial break at atom {i+1}"

    # Validate 3: Residue Identity
    # C3 is the first atom of the patch; verify it was rebranded
    c3_atom = next(at for at in stitched_mol.atoms if at.atom == "C3")
    assert c3_atom.res == "LYS" and c3_atom.res_n == 188, "Patch rebranding failed."

    print("="*60)
    print("RESULT: ALL TOPOLOGY TESTS PASSED")
    print("="*60 + "\n")

if __name__ == "__main__":
    test_full_system_integration()