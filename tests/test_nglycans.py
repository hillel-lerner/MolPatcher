import os
import sys
import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from mol_patcher.pdb_io import PdbParser, PdbBuilder
from mol_patcher.geometry import PatchAligner
from mol_patcher.mol_record import Mol
from mol_patcher.geometry import MolGraph
from mol_patcher.utilities import get_distance


if __name__ == "__main__":
    print("Testing Glycan Bond Alignment...\n")

    pdb_dir = os.path.join(project_root, 'pdbs')
    patch_pdb = os.path.join(pdb_dir, "glycan_only.pdb")
    base_pdb = os.path.join(pdb_dir, 'no-glyc.pdb')
    out_pdb = os.path.join(pdb_dir, 'aligned_glycan_test.pdb')

    target_res = 67    
    target_chain = "A"  
    bond_length = 1.43  # Standard C-N bond length

    patch_res = "BGLC"
    patch_res_num = 1

    print(f"Loading {base_pdb} and {patch_pdb}...")

    parser = PdbParser()
    _, base_records, _ = parser.read_file(base_pdb)
    _, patch_records, _ = parser.read_file(patch_pdb)

    base_records = parser.fix_chain_id(base_records)

    # Protein: Base Parent (CG) and Base Anchor (ND2)
    at1_base_parent = next(r for r in base_records if r.res_seq == target_res and r.chain.strip() == target_chain and r.name.strip() == "CG")
    at2_base_anchor = next(r for r in base_records if r.res_seq == target_res and r.chain.strip() == target_chain and r.name.strip() == "ND2")
    
    # Glycan: Patch Anchor (C1) and Patch Leaving Group (OH leaving group)
    at3_patch_anchor = next(r for r in patch_records if r.name.strip() == "C1")
    at4_patch_leaving = next(r for r in patch_records if r.name.strip() == "O1") 
    at5_patch_leaving = next(r for r in patch_records if r.name.strip() == "HO1") 

    print(f"Aligning glycan {at3_patch_anchor.name} to protein {at2_base_anchor.name}...")

    aligner = PatchAligner(patch_records, [], [])

    aligned_patch = aligner.align_single_bond(
        atoms=patch_records,
        at1=at1_base_parent,
        at2=at2_base_anchor,
        at3=at3_patch_anchor,
        at4=at4_patch_leaving,
        target_bond_length=bond_length
    )

    # Dihedral Test
    try:
        p1_base_parent = next(r for r in base_records if r.res_seq == target_res and r.chain.strip() == target_chain and r.name.strip() == "CB")
        p4_patch_child = next(r for r in patch_records if r.name.strip() == "C2")
        aligned_patch = aligner.set_junction_dihedral(aligned_patch, p1_base_parent, at2_base_anchor, at3_patch_anchor, p4_patch_child, 180)
        print("Dihedral applied successfully.")
    except StopIteration:
        print("Warning: Reference atoms for dihedral not found. Skipping rotation...")

    patch_atoms_to_delete = ["O1", "HO1"]
    filtered_patch = [atom for atom in aligned_patch if atom.name.strip() not in patch_atoms_to_delete]

    base_atoms_to_delete = ["HD22"]
    filtered_base = [atom for atom in base_records if atom.name.strip() not in base_atoms_to_delete]

    combined_records = filtered_base + filtered_patch

    builder = PdbBuilder(out_pdb, combined_records)
    builder.write_pdb()
    print(f"Success!")

