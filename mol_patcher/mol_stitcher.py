import os
from dataclasses import replace
from .utilities import get_distance
from .mol_record import Mol, ItpBond, ItpAngle, ItpDih, ItpPair

def get_pfp_pdb():
    from .pdb_io import PdbParser
    cdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    pfp_path = os.path.join(cdir, 'pdbs', 'pfp_patch_new.pdb')
    _, pfp_atoms, _, _ = PdbParser.read_file(pfp_path)
    return pfp_atoms

def get_pfp_itp(mol_obj):
    cdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    itp_path = os.path.join(cdir, 'itps', 'pfp_patch_new.itp')
    mol_obj.load_itp(itp_path)

def delete_atoms(base_records, target_anchors, extra_deletions=None):
    """Identifies coordinate atoms to delete while protecting anchors."""
    if extra_deletions is None:
        extra_deletions = []
    
    atoms_to_delete = list(extra_deletions)
    distXH = 1.15  
    protected_names = [a.name for a in target_anchors]

    for atom in base_records:
        if atom.name in protected_names:
            continue
        if atom.name.startswith('H'):
            for anchor in target_anchors:
                dist = get_distance([atom.x, atom.y, atom.z], [anchor.x, anchor.y, anchor.z])
                if dist <= distXH:
                    if atom not in atoms_to_delete:
                        atoms_to_delete.append(atom)
                    break 
    return atoms_to_delete

def stitch_molecules(
        base_mol, 
        aligned_patch_atoms, 
        target_reference, 
        target_anchors, 
        patch_mol, 
        patch_anchor_names=["N", "C10", "C11"],
        patch_bridge_name="C7"):
    
    anchor_serial = target_anchors[2].serial
    anchor_res_seq = target_anchors[2].res_seq
    anchor_name = target_anchors[2].name.strip()

    # Handle Deletions
    protein_h_deletions = delete_atoms(base_mol.records, [target_anchors[2]])
    patch_overlap_anchors = [next(a for a in patch_mol.records if a.name.strip() == name) for name in patch_anchor_names]
    patch_to_delete = delete_atoms(patch_mol.records, patch_overlap_anchors)

    # Filter Lists
    final_records = [r for r in base_mol.records if r not in protein_h_deletions]
    
    # Disconnect ITP filtering from PDB index. Identify deleted atoms by matching residue number and atom name.
    deleted_identifiers = [(r.res_seq, r.name.strip()) for r in protein_h_deletions]
    final_atoms = [a for a in base_mol.atoms if (a.res_n, a.atom.strip()) not in deleted_identifiers]
    
    filter_patch_records, filter_patch_atoms = [], []
    dynamic_seg_id = target_anchors[0].seg_id

    for i, record in enumerate(aligned_patch_atoms):
        if record not in patch_to_delete and record not in patch_overlap_anchors:
            filter_patch_records.append(replace(record, res_name=target_reference.res_name, 
                                    chain=target_reference.chain, res_seq=target_reference.res_seq, 
                                    seg_id=dynamic_seg_id))
            filter_patch_atoms.append(replace(patch_mol.atoms[i], res_n=target_reference.res_seq, 
                                    res=target_reference.res_name))

    # Apply Offset to Patch
    offset = 100000
    for atom in filter_patch_atoms:
        atom.number += offset

    patch_bonds = [replace(b, a1=b.a1+offset, a2=b.a2+offset) for b in patch_mol.bonds]
    patch_pairs = [replace(p, a1=p.a1+offset, a2=p.a2+offset) for p in patch_mol.pairs]
    patch_angles = [replace(ang, a1=ang.a1+offset, a2=ang.a2+offset, a3=ang.a3+offset) for ang in patch_mol.angles]
    patch_dihs = [replace(d, a1=d.a1+offset, a2=d.a2+offset, a3=d.a3+offset, a4=d.a4+offset) for d in patch_mol.dihs]

    # Splice Atoms Independently
    insert_idx_records = next(i for i, r in enumerate(final_records) if r.serial == anchor_serial) + 1
    final_records[insert_idx_records:insert_idx_records] = filter_patch_records
    
    # Calculate insertion point for the ITP list based on the target residue and atom name.
    insert_idx_atoms = next(i for i, a in enumerate(final_atoms) if a.res_n == anchor_res_seq and a.atom.strip() == anchor_name) + 1
    final_atoms[insert_idx_atoms:insert_idx_atoms] = filter_patch_atoms

    # Build Initial Object
    stitched_mol = Mol(base_mol.name, final_records, final_atoms, 
                    base_mol.bonds + patch_bonds, base_mol.pairs + patch_pairs, 
                    base_mol.angles + patch_angles, base_mol.dihs + patch_dihs)

    # Add Junction Interactions
    patch_bridge_idx = next(a for a in patch_mol.atoms if a.atom.strip() == patch_bridge_name.strip()).number + offset
    
    # Fetch native ITP numbers
    anch_0_itp = next(a.number for a in base_mol.atoms if a.res_n == target_anchors[0].res_seq and a.atom.strip() == target_anchors[0].name.strip())
    anch_1_itp = next(a.number for a in base_mol.atoms if a.res_n == target_anchors[1].res_seq and a.atom.strip() == target_anchors[1].name.strip())
    anch_2_itp = next(a.number for a in base_mol.atoms if a.res_n == anchor_res_seq and a.atom.strip() == anchor_name)

    stitched_mol.bonds.append(ItpBond(anch_2_itp, patch_bridge_idx, 1))
    stitched_mol.angles.append(ItpAngle(anch_0_itp, anch_2_itp, patch_bridge_idx, 5))
    stitched_mol.dihs.append(ItpDih(anch_1_itp, anch_0_itp, anch_2_itp, patch_bridge_idx, 9))
    stitched_mol.pairs.append(ItpPair(anch_1_itp, patch_bridge_idx, 1))

    # Reindex and Sort
    stitched_mol.reindex()

    stitched_mol.bonds.sort(key=lambda x: x.a1)
    stitched_mol.pairs.sort(key=lambda x: x.a1)
    stitched_mol.angles.sort(key=lambda x: x.a1)
    stitched_mol.dihs.sort(key=lambda x: x.a1)

    return stitched_mol
