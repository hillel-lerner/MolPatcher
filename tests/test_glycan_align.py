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
    print("Testing Glycan Alignment...\n")

    pdb_dir = os.path.join(project_root, 'pdbs')
    proa = os.path.join(pdb_dir, "no-glyc.pdb")
    glyc = os.path.join(pdb_dir, 'm9.pdb')
    out_pdb = os.path.join(pdb_dir, 'aligned_test.pdb')


    glyc_res = 'BGL'
    proa_res = 'ASN'
    proa_chain = 'A'
    # glyc_chain =  

    glyc_res_no = 1
    proa_res_no = 67
    bond_length = 1.43  # Standard C-N bond length
    

    parser = PdbParser()
    _, proa_records, _ = parser.read_file(proa)
    _, glyc_records, _ = parser.read_file(glyc)

    at1_pep_parent = next(r for r in proa_records if r.name.strip() == "CG" and r.res_name == proa_res and r.res_seq == proa_res_no)
    at2_pep_anchor = next(r for r in proa_records if r.name.strip() == "ND2" and r.res_name == proa_res and r.res_seq == proa_res_no)

    # Glycan: Patch Anchor (C1) and Patch Leaving Group (O1)
    at3_glyc_anchor = next(r for r in glyc_records if r.name.strip() == "C1" and r.res_name == glyc_res and r.res_seq == glyc_res_no)
    at4_glyc_leaving = next(r for r in glyc_records if r.name.strip() == "O1" and r.res_name == glyc_res and r.res_seq == glyc_res_no)

    at_od1 = next(r for r in proa_records if r.name.strip() == "OD1" and r.res_name == proa_res and r.res_seq == proa_res_no)


    # repulsion coordinates
    base_coords = np.array([[r.x, r.y, r.z] for r in proa_records])

    # Calculate mean across the rows (axis=0) to get the centroid
    base_centroid = np.mean(base_coords, axis=0)

    aligner = PatchAligner(glyc_records, [], [])

    aligned_glyc = aligner.align_single_bond(
        atoms=glyc_records,
        at1=at1_pep_parent,
        at2=at2_pep_anchor,
        at3=at3_glyc_anchor,
        at4=at4_glyc_leaving,
        target_bond_length=bond_length
    )

    target_omega = 180.0
    
    aligned_glyc = aligner.set_junction_dihedral(
        atoms=aligned_glyc, 
        p1=at_od1, 
        p2=at1_pep_parent, # This is your CG
        p3=at2_pep_anchor, # This is your ND2
        p4=at3_glyc_anchor, # This is your C1
        target_dih=target_omega
    )
    
    print(f"Success: Omega_N locked to {target_omega}° (trans-planar configuration).")

    p_ref = next(r for r in proa_records if r.name.strip() == "CB" and r.res_name == proa_res and r.res_seq == proa_res_no)

    # Apply the angle adjustment
    aligned_glyc = aligner.set_junction_angle(
        atoms=aligned_glyc,
        p1=at1_pep_parent,   # ASN CG
        p2=at2_pep_anchor,   # ASN ND2
        p3=at3_glyc_anchor,  # Glycan C1
        p_ref=p_ref,         # ASN CB (defines the rotation plane)
        target_angle=120.0,  # 120.0 for sp2 Amide Nitrogen
        repulsion_coord=base_centroid 
    )

    

    p4_glyc_child = next(r for r in glyc_records if r.name.strip() == "O5")
    target_angle = -120 

    aligned_glyc = aligner.set_junction_dihedral(
        atoms=aligned_glyc, 
        p1=at1_pep_parent, 
        p2=at2_pep_anchor, 
        p3=at3_glyc_anchor, 
        p4=p4_glyc_child, 
        target_dih=target_angle
    )
    print(f"Dihedral twisted to {target_angle} degrees.")

    proa_atoms_to_delete = ["HD21"] 
    filtered_proa = [atom for atom in proa_records if atom.name.strip() not in proa_atoms_to_delete]

    glyc_atoms_to_delete = ["O1", "HO1"]
    filtered_glyc = [atom for atom in aligned_glyc if atom.name.strip() not in glyc_atoms_to_delete]

    combined_records = filtered_proa + filtered_glyc
    
    for i, atom in enumerate(combined_records):
        atom.serial = i + 1

    builder = PdbBuilder(out_pdb, combined_records)
    builder.write_pdb()
    print(f"Success!")
