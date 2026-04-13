from typing import List, Optional
from dataclasses import dataclass, field, replace
from .mol_record import Mol

def reindex_topology(mol, remap_dict=None):
    """
    Synchronizes atom numbers and all bonded interactions (bonds, pairs, angles, dihedrals)
    to ensure sequential numbering after molecular surgery.
    
    :param mol: (Mol) The molecule object to be reindexed.
    :param remap_dict: (dict, optional) A dictionary mapping deleted or merged atom indices 
                    to their new preserved index (e.g., mapping a patch nitrogen to a protein nitrogen).
    :return: (dict) index_map containing the translation from old indices to new sequential indices.
    """
    
    index_map = {}
    remap_dict = remap_dict or {}

    for i, atom in enumerate(mol.atoms):
        new_nr = i + 1
        index_map[atom.number] = new_nr
        atom.number = new_nr
        atom.cgnr = new_nr  
        
    for i, record in enumerate(mol.records):
        record.serial = i + 1

    # Update Bonds
    new_bonds = []
    seen_bonds = set()
    for bond in mol.bonds:
        # Check if either atom was merged/deleted and needs to be pointed to a new anchor
        target_a1 = remap_dict.get(bond.a1, bond.a1)
        target_a2 = remap_dict.get(bond.a2, bond.a2)
        try:
            # Translate the old or remapped index into the new sequential index
            new_a1 = index_map[target_a1]
            new_a2 = index_map[target_a2]
            
            # Use sorted tuples to prevent duplicate interactions backwards and forwards
            signature = tuple(sorted([new_a1, new_a2]))
            if signature not in seen_bonds:
                seen_bonds.add(signature)
                new_bonds.append(replace(bond, a1=new_a1, a2=new_a2))
        except KeyError:
            continue
    mol.bonds = new_bonds

    # Update Pairs
    new_pairs = []
    seen_pairs = set()
    for pair in mol.pairs:
        target_a1 = remap_dict.get(pair.a1, pair.a1)
        target_a2 = remap_dict.get(pair.a2, pair.a2)
        try:
            new_a1 = index_map[target_a1]
            new_a2 = index_map[target_a2]
            
            signature = tuple(sorted([new_a1, new_a2]))
            if signature not in seen_pairs:
                seen_pairs.add(signature)
                new_pairs.append(replace(pair, a1=new_a1, a2=new_a2))
        except KeyError:
            continue
    mol.pairs = new_pairs

    # Update Angles
    new_angles = []
    seen_angles = set()
    for ang in mol.angles:
        target_a1 = remap_dict.get(ang.a1, ang.a1)
        target_a2 = remap_dict.get(ang.a2, ang.a2)
        target_a3 = remap_dict.get(ang.a3, ang.a3)
        try:
            new_a1 = index_map[target_a1]
            new_a2 = index_map[target_a2]
            new_a3 = index_map[target_a3]
            
            # Check forwards and backwards
            sig1 = (new_a1, new_a2, new_a3, ang.type)
            sig2 = (new_a3, new_a2, new_a1, ang.type)
            
            if sig1 not in seen_angles and sig2 not in seen_angles:
                seen_angles.add(sig1)
                new_angles.append(replace(ang, a1=new_a1, a2=new_a2, a3=new_a3))
        except KeyError:
            continue
    mol.angles = new_angles

    # Update Dihedrals
    new_dihs = []
    seen_dihs = set()
    for dih in mol.dihs:
        target_a1 = remap_dict.get(dih.a1, dih.a1)
        target_a2 = remap_dict.get(dih.a2, dih.a2)
        target_a3 = remap_dict.get(dih.a3, dih.a3)
        target_a4 = remap_dict.get(dih.a4, dih.a4)
        try:
            new_a1 = index_map[target_a1]
            new_a2 = index_map[target_a2]
            new_a3 = index_map[target_a3]
            new_a4 = index_map[target_a4]
            
            # Check forwards and backwards
            sig1 = (new_a1, new_a2, new_a3, new_a4, dih.type)
            sig2 = (new_a4, new_a3, new_a2, new_a1, dih.type)
            
            if sig1 not in seen_dihs and sig2 not in seen_dihs:
                seen_dihs.add(sig1)
                new_dihs.append(replace(dih, a1=new_a1, a2=new_a2, a3=new_a3, a4=new_a4))
        except KeyError:
            continue
    mol.dihs = new_dihs

    return index_map