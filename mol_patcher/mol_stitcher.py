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
    
    """
    Inserts patch atoms immediately after the anchor NZ atom
    Calls Mol.reindex() to update the atom numbers for everything at the insertion point and onwards.
    Modifies the subsequent calls to each of the updated atoms for the itp files.
    """
    anchor_serial = target_anchors[2].serial  # anchor serial is NZ on the protein

    # Filter existing atoms
    protein_h_deletions = delete_atoms(base_mol.records, [target_anchors[2]])
    
    patch_overlap_anchors = [
        next(a for a in patch_mol.records if a.name.strip() == name) 
        for name in patch_anchor_names
    ]

    patch_to_delete = delete_atoms(patch_mol.records, patch_overlap_anchors)

    # Filter protein atoms and find the exact insertion index
    final_records = [r for r in base_mol.records if r not in protein_h_deletions]
    final_atoms = [a for i, a in enumerate(base_mol.atoms) if base_mol.records[i] not in protein_h_deletions]
    
    insert_idx = next(i for i, r in enumerate(final_records) if r.serial == anchor_serial) + 1

    # Filter and re-tag the patch atoms
    dynamic_seg_id = target_anchors[0].seg_id
    filter_patch_records, filter_patch_atoms = [], []

    for i, record in enumerate(aligned_patch_atoms):
        if record not in patch_to_delete and record not in patch_overlap_anchors:
            filter_patch_records.append(replace(record, res_name=target_reference.res_name, 
                                    chain=target_reference.chain, res_seq=target_reference.res_seq, 
                                    seg_id=dynamic_seg_id))
            filter_patch_atoms.append(replace(patch_mol.atoms[i], res_n=target_reference.res_seq, 
                                    res=target_reference.res_name))

    # Splice patch into the middle of the lists
    final_records[insert_idx:insert_idx] = filter_patch_records
    final_atoms[insert_idx:insert_idx] = filter_patch_atoms

    # Assemble and add topological interactions at the NZ-C7 junction)
    stitched_mol = Mol(base_mol.name, final_records, final_atoms, 
                    base_mol.bonds + patch_mol.bonds, base_mol.pairs + patch_mol.pairs, 
                    base_mol.angles + patch_mol.angles, base_mol.dihs + patch_mol.dihs)

    patch_bridge_idx = next(a for a in patch_mol.atoms if a.atom.strip() == patch_bridge_name.strip()).number

    stitched_mol.bonds.append(ItpBond(anchor_serial, patch_bridge_idx, 1))
    stitched_mol.angles.append(ItpAngle(target_anchors[0].serial, anchor_serial, patch_bridge_idx, 5))
    stitched_mol.dihs.append(ItpDih(target_anchors[1].serial, target_anchors[0].serial, anchor_serial, patch_bridge_idx, 9))
    stitched_mol.pairs.append(ItpPair(target_anchors[1].serial, patch_bridge_idx, 1))

    # Synchronize the pdb + itp information
    stitched_mol.reindex()
    return stitched_mol
