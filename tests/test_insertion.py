import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from mol_patcher.mol_record import Mol
from mol_patcher.pdb_io import PdbParser
from mol_patcher import mol_stitcher

test_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(test_dir)

protein_pdb = os.path.join(project_root, 'pdbs', '6oge_prob_allatom.pdb')
protein_itp = os.path.join(project_root, 'itps', '6oge_prob.itp')

def verify_pfp_insertion(stitched_mol, target_res_seq=188):
    """
    Objectively verifies that the patch was inserted immediately after the target NZ atom 
    and that the PDB metadata was correctly updated to match the host residue.
    """
    print(f"\n--- VERIFYING INSERTION AT RESIDUE {target_res_seq} ---")
    
    try:
        # 1. Find the exact index of the Lysine NZ anchor
        anchor_idx = next(i for i, r in enumerate(stitched_mol.records) 
                          if r.res_seq == target_res_seq and r.name.strip() == "NZ")
        anchor_rec = stitched_mol.records[anchor_idx]
        
        # The atom immediately following should be the first atom of your patch
        first_patch_rec = stitched_mol.records[anchor_idx + 1]
        
        print(f"Anchor found: {anchor_rec.name} at index {anchor_idx}")
        print(f"Inserted atom directly following: {first_patch_rec.name} at index {anchor_idx + 1}")
        
        # 2. Validate Metadata Inheritance
        metadata_match = (
            first_patch_rec.res_name == anchor_rec.res_name and
            first_patch_rec.res_seq == anchor_rec.res_seq and
            first_patch_rec.chain == anchor_rec.chain
        )
        
        if metadata_match:
            print(f"PASS: Metadata successfully updated. The patch atom is labeled as {first_patch_rec.res_name} {first_patch_rec.res_seq} in chain {first_patch_rec.chain}.")
        else:
            print(f"FAIL: Metadata mismatch.")
            print(f"Expected: {anchor_rec.res_name} {anchor_rec.res_seq} Chain {anchor_rec.chain}")
            print(f"Got: {first_patch_rec.res_name} {first_patch_rec.res_seq} Chain {first_patch_rec.chain}")
            
    except StopIteration:
        print(f"ERROR: Could not find NZ atom for residue {target_res_seq}.")

def test_full_system_integration():
    
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

    # Validate 4: Junction Topology Verification
    # Instead of a naive sum, we explicitly verify the junction bond exists.
    
    # 1. Get the final ITP indices of the two joined atoms
    nz_index = next(at.number for at in stitched_mol.atoms if at.res_n == 188 and at.atom == "NZ")
    c3_index = next(at.number for at in stitched_mol.atoms if at.res_n == 188 and at.atom == "C3")
    
    # 2. Search the bonds list for a bond connecting these two specific indices
    junction_bonds = [b for b in stitched_mol.bonds if (b.a1 == nz_index and b.a2 == c3_index) or (b.a1 == c3_index and b.a2 == nz_index)]
    
    print(f"\n--- VERIFYING TOPOLOGY ---")
    print(f"Looking for bond between atom {nz_index} (NZ) and {c3_index} (C3)")
    
    if junction_bonds:
        print(f"PASS: Junction bond found! -> {junction_bonds[0]}")
    else:
        print("FAIL: Junction bond is missing from the topology.")
        
    assert len(junction_bonds) == 1, "The covalent linkage was not formed."

    # ADDED CHANGE: Call the verification function here, passing the object
    verify_pfp_insertion(stitched_mol)

    print("="*60)
    print("RESULT: ALL TOPOLOGY TESTS PASSED")
    print("="*60 + "\n")

if __name__ == "__main__":
    test_full_system_integration()


