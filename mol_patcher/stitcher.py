import os
import numpy as np
from dataclasses import replace
from .geometry import MolGraph
from .utilities import get_distance
from .topology_tools import reindex_topology
from .mol_record import Mol, ItpBond, ItpAngle, ItpDih, ItpPair, ItpAtom
from .pdb_io import PdbParser


class Patchloader:
    """Loads the ligand patch pdb and itp files into mol objects for processing"""

    def __init__(self):
        self.f = os.path.abspath(__file__)
        self.mol_patcher_dir = os.path.dirname(self.f)
        self.root = os.path.dirname(self.mol_patcher_dir)
        self.pdb_dir = os.path.join(self.root, "pdbs")
        self.top_dir = os.path.join(self.root, "topology_files")
    
    def get_pfp_pdb(self):
        pfp_path = os.path.join(self.pdb_dir, "pfp_patch_new.pdb")
        _, pfp_atoms, _, _ = PdbParser.read_file(pfp_path)
        return pfp_atoms

    
    def get_pfp_itp(self, mol_obj):
        # Loads the ITP topology into an existing Mol object.
        itp_path = os.path.join(self.top_dir, "pfp_patch_new.itp")
        mol_obj.load_itp(itp_path)


class Stitcher:
    """Performs the molecular surgury to attach the base molecule and the patch molecule"""
    def __init__(self, base_mol, patch_mol, target_res_id):
        self.base = base_mol
        self.patch = patch_mol
        self.res_id = target_res_id

        # State tracking
        self.offset = 100000
        self.base_deletions = []
        self.patch_deletions = []
        self.junction_interactions = []

    def get_bonded_hydrogens(self, mol, mol_graph, anchors):
        """Uses the pre-built molecular graph to instantly find hydrogens bonded to the anchor atoms."""
        h_to_delete = []
        for anchor in anchors:
            # Find the exact index of the anchor in the molecule's records using its serial number
            anchor_idx = next(i for i, r in enumerate(mol.records) if r.serial == anchor.serial)
            
            # Fetch all bonded neighbors from the networkx graph
            for n_idx in mol_graph.nx_graph.neighbors(anchor_idx):
                neighbor_record = mol.records[n_idx]
                if neighbor_record.name.strip().startswith('H'):
                    if neighbor_record not in h_to_delete:
                        h_to_delete.append(neighbor_record)
        return h_to_delete
    

    def build_bridge_topol(self, stitched_mol, stitched_graph, nz_pdb_idx, c7_pdb_idx):
        """
        Uses the connectivity matrix to check for missing topology data in the patched molecule
        """
        nz_neighbors = list(stitched_graph.nx_graph.neighbors(nz_pdb_idx))
        c7_neighbors = list(stitched_graph.nx_graph.neighbors(c7_pdb_idx))
        nz_exclusive = [n for n in nz_neighbors if n != c7_pdb_idx]
        c7_exclusive = [n for n in c7_neighbors if n != nz_pdb_idx]

        # Translate a PDB graph index into an ITP number
        def to_itp(pdb_idx):
            record = stitched_mol.records[pdb_idx]
            return next(a.number for a in stitched_mol.atoms if a.res_n == record.res_seq and a.atom.strip() == record.name.strip())

        # Translate the bridge atom indices
        nz_itp = to_itp(nz_pdb_idx)
        c7_itp = to_itp(c7_pdb_idx)

        angles_to_check = []
        for n in nz_exclusive:
            angles_to_check.append((to_itp(n), nz_itp, c7_itp))
        for c in c7_exclusive:
            angles_to_check.append((nz_itp, c7_itp, to_itp(c)))

        dihedrals_to_check = []
        for n in nz_exclusive:
            for c in c7_exclusive:
                dihedrals_to_check.append((to_itp(n), nz_itp, c7_itp, to_itp(c)))

        print("Angles across junction (ITP Numbers):", angles_to_check)
        print("Dihedrals across junction (ITP Numbers):", dihedrals_to_check)

        return angles_to_check, dihedrals_to_check
    

    def stitch_molecules(self, aligned_patch_atoms, target_reference, target_anchors, patch_anchor_names=["N", "C10", "C11"], patch_bridge_name="C7"):
        anchor_serial = target_anchors[2].serial
        anchor_res_seq = target_anchors[2].res_seq
        anchor_name = target_anchors[2].name.strip()

        print("Building molecular graphs for deletion analysis...")
        base_graph = MolGraph(self.base)
        patch_graph = MolGraph(self.patch)

        self.base_deletions = self.get_bonded_hydrogens(self.base, base_graph, [target_anchors[2]])

        hz_candidates = [a for a in self.base_deletions if a.name.strip() in ['HZ1', 'HZ2', 'HZ3']]
        if len(hz_candidates) > 0:
            kept_h = hz_candidates[0]
            self.base_deletions.remove(kept_h)
            print(f"Preserved {kept_h.name} for the junction.")

            
        patch_overlap_anchors = [next(a for a in self.patch.records if a.name.strip() == name) for name in patch_anchor_names]
        self.patch_deletions = self.get_bonded_hydrogens(self.patch, patch_graph, patch_overlap_anchors)

        print(f"Final kill-list: {len(self.base_deletions)} base atoms, {len(self.patch_deletions)} patch atoms.")

        # Filter Lists
        final_records = [r for r in self.base.records if r not in self.base_deletions]
        
        # Disconnect ITP filtering from PDB index. Identify deleted atoms by matching residue number and atom name.
        deleted_identifiers = [(r.res_seq, r.name.strip()) for r in self.base_deletions]
        final_atoms = [a for a in self.base.atoms if (a.res_n, a.atom.strip()) not in deleted_identifiers]
        
        filter_patch_records, filter_patch_atoms = [], []
        dynamic_seg_id = target_anchors[0].seg_id

        for i, record in enumerate(aligned_patch_atoms):
            if record not in self.patch_deletions and record not in patch_overlap_anchors:
                filter_patch_records.append(replace(record, res_name=target_reference.res_name, 
                                        chain=target_reference.chain, res_seq=target_reference.res_seq, 
                                        seg_id=dynamic_seg_id))
                filter_patch_atoms.append(replace(self.patch.atoms[i], res_n=target_reference.res_seq, 
                                        res=target_reference.res_name))

        for atom in filter_patch_atoms:
            atom.number += self.offset

        patch_bonds = [replace(b, a1=b.a1+self.offset, a2=b.a2+self.offset) for b in self.patch.bonds]
        patch_pairs = [replace(p, a1=p.a1+self.offset, a2=p.a2+self.offset) for p in self.patch.pairs]
        patch_angles = [replace(ang, a1=ang.a1+self.offset, a2=ang.a2+self.offset, a3=ang.a3+self.offset) for ang in self.patch.angles]
        patch_dihs = [replace(d, a1=d.a1+self.offset, a2=d.a2+self.offset, a3=d.a3+self.offset, a4=d.a4+self.offset) for d in self.patch.dihs]

        # Splice Atoms Independently
        insert_idx_records = next(i for i, r in enumerate(final_records) if r.serial == anchor_serial) + 1
        final_records[insert_idx_records:insert_idx_records] = filter_patch_records
        
        # Calculate insertion point for the ITP list based on the target residue and atom name.
        insert_idx_atoms = next(i for i, a in enumerate(final_atoms) if a.res_n == anchor_res_seq and a.atom.strip() == anchor_name) + 1
        final_atoms[insert_idx_atoms:insert_idx_atoms] = filter_patch_atoms

        # Build Initial Object
        stitched_mol = Mol(self.base.name, final_records, final_atoms, 
                        self.base.bonds + patch_bonds, self.base.pairs + patch_pairs, 
                        self.base.angles + patch_angles, self.base.dihs + patch_dihs)

        # Add Junction Interactions (1-2 Bond Only)
        patch_bridge_idx = next(a for a in self.patch.atoms if a.atom.strip() == patch_bridge_name.strip()).number + self.offset
        anch_2_itp = next(a.number for a in self.base.atoms if a.res_n == anchor_res_seq and a.atom.strip() == anchor_name)
        
        stitched_mol.bonds.append(ItpBond(anch_2_itp, patch_bridge_idx, 1))

        # Remap topologies: Map N -> NZ 
        patch_n_idx = next(a.number for a in self.patch.atoms if a.atom.strip() == "N") + self.offset
        remap_dict = {
            patch_n_idx: anch_2_itp   
        }
        
        # Reindex and Sort
        reindex_topology(stitched_mol, remap_dict)
        
        # Build graph on the final, cleanly-numbered molecule
        stitched_graph = MolGraph(stitched_mol)

        # Find the PDB graph indices of the NZ and C7 atoms
        nz_idx = next(i for i, r in enumerate(stitched_mol.records) if r.res_seq == anchor_res_seq and r.name.strip() == anchor_name)
        c7_idx = next(i for i, r in enumerate(stitched_mol.records) if r.res_seq == target_reference.res_seq and r.name.strip() == patch_bridge_name.strip())

        # Call the checker to get the missing parameters
        angles_to_check, dihedrals_to_check = self.build_bridge_topol(stitched_mol, stitched_graph, nz_idx, c7_idx)
        
        # --- AUTOMATED APPENDS ---
        # Unpack 3-item angle tuples
        for a1, a2, a3 in angles_to_check:
            stitched_mol.angles.append(ItpAngle(a1, a2, a3, 5))

        # Unpack the 4-item dihedral tuples
        for a1, a2, a3, a4 in dihedrals_to_check:
            stitched_mol.dihs.append(ItpDih(a1, a2, a3, a4, 9))
            stitched_mol.pairs.append(ItpPair(a1, a4, 1)) # 1-4 pairs always use the outer atoms
            
        # Sort Final Lists
        stitched_mol.bonds.sort(key=lambda x: x.a1)
        stitched_mol.pairs.sort(key=lambda x: x.a1)
        stitched_mol.angles.sort(key=lambda x: x.a1)
        stitched_mol.dihs.sort(key=lambda x: x.a1)

        type_map={
            "NZ": "NG311",
            "HZ1": "HGPAM1", 
            "HZ2": "HGPAM1", 
            "HZ3": "HGPAM1"
        }  
        
        stitched_mol = self.update_atom_types(stitched_mol, type_map)
        stitched_mol = self.balance_electrostatics(stitched_mol)

        return stitched_mol
    
    def update_atom_types(self, stitched_mol, type_map):
        """Updates specific atom types at the junction (e.g., NZ -> NG311)."""
        for atom in stitched_mol.atoms:
            if atom.res_n == self.res_id and atom.atom.strip() in type_map:
                atom.type = type_map[atom.atom.strip()]
        return stitched_mol
    
    def balance_electrostatics(self, stitched_mol):

        # Update the anchor charges first
        # Dictionary maps the base atom name to its equivalent patch atom name
        anchor_map = {'NZ': 'N', 'CE': 'C10', 'CD': 'C11'}
        
        for base_atom in stitched_mol.atoms:
            if base_atom.res_n == self.res_id and base_atom.atom.strip() in anchor_map:
                patch_equiv_name = anchor_map[base_atom.atom.strip()]
                
                # Look up the new charge from the raw patch molecule
                patch_equiv_atom = next(a for a in self.patch.atoms if a.atom.strip() == patch_equiv_name)
                
                old_charge = float(base_atom.charge)
                base_atom.charge = float(patch_equiv_atom.charge)
                base_atom.comment = f"; old charge: {old_charge:.4f} (anchor update)"


        # Calculate the total global deficit
        current_sum = sum(float(a.charge) for a in stitched_mol.atoms)
        target_sum = round(current_sum)
        delta_q = target_sum - current_sum


        # Distribute the deficit across the sinks to achieve perfect integer
        sinks = ['HD1' , 'HD2', 'HG1', 'HG2', 'H2', 'H3', 'H4', 'H5', 'H6']
        q_per_atom = delta_q / len(sinks)

        for i in stitched_mol.atoms:
            if i.res_n == self.res_id and i.atom.strip() in sinks:
                current_charge = float(i.charge)
                i.charge = current_charge + q_per_atom
                i.comment = f"; old charge: {current_charge:.4f}, delta q = {q_per_atom:.4f}"

        # Final verification print
        total_stitched_charge = sum(float(a.charge) for a in stitched_mol.atoms)
        print(f"Total Stitched Charge: {total_stitched_charge:.6f}")

        return stitched_mol