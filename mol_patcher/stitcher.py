import os
import numpy as np
from dataclasses import replace
from .geometry import MolGraph
from .utilities import get_distance
from .topology_tools import reindex_topology
from .mol_record import Mol, ItpBond, ItpAngle, ItpDih, ItpPair, ItpAtom
from .pdb_io import PdbParser


class Patchloader:
    """
    Locates and loads the standard patch ligand PDB and ITP files into Mol objects for processing.
    """

    def __init__(self):
        self.f = os.path.abspath(__file__)
        self.mol_patcher_dir = os.path.dirname(self.f)
        self.root = os.path.dirname(self.mol_patcher_dir)
        self.pdb_dir = os.path.join(self.root, "pdbs")
        self.top_dir = os.path.join(self.root, "topology_files")
    
    def get_pfp_pdb(self):
        pfp_path = os.path.join(self.pdb_dir, "pfp_patch_new.pdb")
        parser = PdbParser()
        _, pfp_atoms, _ = parser.read_file(pfp_path)
        return pfp_atoms

    
    def get_pfp_itp(self, mol_obj):
        # Loads the ITP topology into an existing Mol object.
        itp_path = os.path.join(self.top_dir, "pfp_patch_new.itp")
        mol_obj.load_itp(itp_path)


class Stitcher:
    """
    Performs the topological and geometric surgery required to attach a patch molecule 
    to a base protein molecule.
    """
    def __init__(self, base_mol, patch_mol, target_res_id, itp_chains=None):
        self.base = base_mol
        self.patch = patch_mol
        self.res_id = target_res_id
        self.itp_chains = itp_chains or []

        # State tracking
        self.offset = 100000
        self.base_deletions = []
        self.patch_deletions = []
        self.junction_interactions = []

    def get_bonded_hydrogens(self, mol, mol_graph, anchors):
        """
        Uses the pre-built molecular connectivity graph to identify hydrogens bonded to specific anchors.
        
        :param mol: (Mol) The molecule object containing the records.
        :param mol_graph: (MolGraph) The connectivity matrix/graph of the molecule.
        :param anchors: (list) PdbRecord objects representing the anchor atoms.
        :return: (list) A list of PdbRecord objects representing the bonded hydrogens to be deleted.
        """
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
    
    def _build_pdb_to_itp_map(self, records, atoms, itp_chains=None):
        """
        Creates a dictionary mapping the Python id() of a PdbRecord
        to its corresponding ItpAtom object. Filters PDB records to 
        match the scope of the provided ITP chains.
        """
        # Build ITP blocks
        itp_blocks = []
        current_itp_res = None
        for a in atoms:
            if a.res_n != current_itp_res:
                itp_blocks.append({})
                current_itp_res = a.res_n
            itp_blocks[-1][a.atom.strip()] = a

        # Build Scoped PDB blocks (based on chains)
        scoped_pdb_blocks = []
        current_pdb_res = None
        
        for r in records:
            # Skip records that aren't part of the ITP file's scope
            if itp_chains and r.chain.strip() not in itp_chains:
                continue
                
            res_id = (r.chain.strip(), r.res_seq, r.res_name.strip())
            if res_id != current_pdb_res:
                scoped_pdb_blocks.append({})
                current_pdb_res = res_id
            scoped_pdb_blocks[-1][r.name.strip()] = r

        # Dynamic sequence alignment
        mapping = {}
        i, j = 0, 0
        
        while i < len(scoped_pdb_blocks) and j < len(itp_blocks):
            p_res = list(scoped_pdb_blocks[i].values())[0].res_name.strip()
            i_res = list(itp_blocks[j].values())[0].res.strip()
            
            # Match first 2 chars (LYS/LYX, HIS/HSD)
            if p_res[:2] == i_res[:2]:
                for atom_name, record_obj in scoped_pdb_blocks[i].items():
                    if atom_name in itp_blocks[j]:
                        mapping[id(record_obj)] = itp_blocks[j][atom_name]
                i += 1
                j += 1
            else:
                # Mismatch if the sequences diverge.
                # To prevent false matches, scan ahead in BOTH files to find a unique fingerprint: 3 consecutive residues that perfectly match.
                found_sync = False
                
                # Check up to 40 residues ahead in both files
                for offset_i in range(40):
                    if i + offset_i >= len(scoped_pdb_blocks): break
                    for offset_j in range(40):
                        if j + offset_j >= len(itp_blocks): break
                        
                        # Validate a 3-residue sequence 
                        sync_length = min(3, len(scoped_pdb_blocks) - (i + offset_i), len(itp_blocks) - (j + offset_j))
                        is_match = True
                        
                        for k in range(sync_length):
                            pk = list(scoped_pdb_blocks[i + offset_i + k].values())[0].res_name.strip()[:2]
                            ik = list(itp_blocks[j + offset_j + k].values())[0].res.strip()[:2]
                            if pk != ik:
                                is_match = False
                                break
                                
                        # If the fingerprint matches, lock the alignment and break
                        if is_match and sync_length > 0:
                            i = i + offset_i
                            j = j + offset_j
                            found_sync = True
                            break
                            
                    if found_sync:
                        break
                        
                # If no sequence fingerprint can be found, safely advance PDB to prevent infinite loops
                if not found_sync:
                    i += 1

        return mapping

    def build_bridge_topol(self, stitched_mol, stitched_graph, nz_pdb_idx, c7_pdb_idx, stitched_map):
        """
        Uses the connectivity matrix to identify and generate missing angle and dihedral 
        topology parameters across the newly formed junction bond.
        
        :param stitched_mol: (Mol) The combined molecule object.
        :param stitched_graph: (MolGraph) The connectivity graph of the combined molecule.
        :param nz_pdb_idx: (int) The Python list index of the base protein anchor (e.g., NZ).
        :param c7_pdb_idx: (int) The Python list index of the patch anchor (e.g., C7).
        :param stitched_map: (dict) The map relating the pdb indices to the itp atom indices
        :return: (tuple) Two lists containing the angle tuples and dihedral tuples to be appended.
        """
        nz_neighbors = list(stitched_graph.nx_graph.neighbors(nz_pdb_idx))
        c7_neighbors = list(stitched_graph.nx_graph.neighbors(c7_pdb_idx))
        nz_exclusive = [n for n in nz_neighbors if n != c7_pdb_idx]
        c7_exclusive = [n for n in c7_neighbors if n != nz_pdb_idx]

        # Translate a PDB graph index into an ITP number
        def to_itp(pdb_idx):
            record = stitched_mol.records[pdb_idx]
            mapped_atom = stitched_map.get(id(record))
            
            if mapped_atom is not None:
                return mapped_atom.number
            
            for a in stitched_mol.atoms:
                if a.res_n == record.res_seq and a.atom.strip() == record.name.strip():
                    return a.number

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

        # print("Angles across junction (ITP Numbers):", angles_to_check)
        # print("Dihedrals across junction (ITP Numbers):", dihedrals_to_check)

        return angles_to_check, dihedrals_to_check
    

    def stitch_molecules(self, aligned_patch_atoms, target_reference, target_anchors, patch_anchor_names=["N", "C10", "C11"], patch_bridge_name="C7"):
        """
        Executes the full attachment pipeline: deletes overlapping atoms, offsets topologies, 
        splices lists, generates junction bonds, and balances electrostatics.
        
        :param aligned_patch_atoms: (list) PdbRecord objects of the patch translated to the target site.
        :param target_reference: (PdbRecord) The primary target anchor used for chain/residue designation.
        :param target_anchors: (list) PdbRecord objects of the protein attachment site.
        :param patch_anchor_names: (list) String names of the patch atoms that will overlap and be deleted.
        :param patch_bridge_name: (str) String name of the patch atom forming the new bond.
        :return: (Mol) The final stitched, reindexed, and charge-balanced molecule.
        """
        
        anchor_serial = target_anchors[2].serial
        anchor_res_seq = target_anchors[2].res_seq
        anchor_name = target_anchors[2].name.strip()

        base_map = self._build_pdb_to_itp_map(self.base.records, self.base.atoms, itp_chains=self.itp_chains)

        old_itp_to_pdb = {}
        for record in self.base.records:
            itp_obj = base_map.get(id(record))
            if itp_obj is not None:
                old_itp_to_pdb[itp_obj.number] = record

        print("Building molecular graphs for deletion analysis...")
        base_graph = MolGraph(self.base, itp_to_pdb_map=old_itp_to_pdb)
        # Patch doesn't have a map, so it falls back to distance checks
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

        # Append patch atoms with the LYX designation and record old residue number
        for i, record in enumerate(aligned_patch_atoms):
            if record not in self.patch_deletions and record not in patch_overlap_anchors:
                filter_patch_records.append(replace(record, res_name="LYX", 
                                        chain=target_reference.chain, res_seq=target_reference.res_seq, 
                                        seg_id=dynamic_seg_id))
                original_atom = self.patch.atoms[i]
                old_res_n = original_atom.res_n
                new_atom = replace(original_atom, res_n=target_reference.res_seq, res="LYX")

                if old_res_n != target_reference.res_seq:
                    existing_comment = new_atom.comment.strip() if hasattr(new_atom, 'comment') and new_atom.comment else ""
                    new_atom.comment = f"{existing_comment} ; old_res: {old_res_n}".strip()
                    
                filter_patch_atoms.append(new_atom)

        # Rename the preserved base protein atoms at the target site to LYX
        for r in final_records:
            if r.res_seq == target_reference.res_seq and r.chain == target_reference.chain and r.res_name.strip() == "LYS":
                r.res_name = "LYX"

        for atom in filter_patch_atoms:
            atom.number += self.offset

        patch_bonds = [replace(b, a1=b.a1+self.offset, a2=b.a2+self.offset) for b in self.patch.bonds]
        patch_pairs = [replace(p, a1=p.a1+self.offset, a2=p.a2+self.offset) for p in self.patch.pairs]
        patch_angles = [replace(ang, a1=ang.a1+self.offset, a2=ang.a2+self.offset, a3=ang.a3+self.offset) for ang in self.patch.angles]
        patch_dihs = [replace(d, a1=d.a1+self.offset, a2=d.a2+self.offset, a3=d.a3+self.offset, a4=d.a4+self.offset) for d in self.patch.dihs]

        # Splice Atoms Independently
        insert_idx_records = next(i for i, r in enumerate(final_records) if r.serial == anchor_serial) + 1
        final_records[insert_idx_records:insert_idx_records] = filter_patch_records
        
        # Sync base ITP residue numbers to the PDB and track old numbers
        for a in final_atoms:
            rec = old_itp_to_pdb.get(a.number)
            if rec:
                new_res_n = rec.res_seq
                
                # Extract original number from clean_itp_file if it exists
                true_old_res = str(a.res_n)
                if 'old_res:' in getattr(a, 'comment', ''):
                    true_old_res = a.comment.split('old_res:')[-1].strip()
                    
                # Sync the actual integer topology number to the PDB
                if a.res_n != new_res_n:
                    a.res_n = new_res_n
                    
                # If the original number perfectly matches the final PDB number, wipe the comment
                if str(new_res_n) == true_old_res:
                    a.comment = ""
                else:
                    # Only print if it changed
                    a.comment = f"; old_res: {true_old_res}"
            
                if rec.res_seq == target_reference.res_seq and rec.chain.strip() == target_reference.chain.strip():
                    if a.res.strip() == "LYS":
                        a.res = "LYX"
        
        anch_2_itp_obj = base_map.get(id(target_anchors[2]))
        if anch_2_itp_obj is None:
            for a in self.base.atoms:
                if a.res_n == target_anchors[2].res_seq and a.atom.strip() == target_anchors[2].name.strip():
                    anch_2_itp_obj = a
                    break

        insert_idx_atoms = final_atoms.index(anch_2_itp_obj) + 1
        
        final_atoms[insert_idx_atoms:insert_idx_atoms] = filter_patch_atoms

        # Build Initial Object
        stitched_mol = Mol(self.base.name, self.base.moltype_section, final_records, final_atoms, 
                        patch_bonds, self.base.pairs + patch_pairs, # <-- removed self.base.bonds
                        self.base.angles + patch_angles, self.base.dihs + patch_dihs)

        # Add Junction Interactions (1-2 Bond Only)
        patch_bridge_idx = next(a for a in self.patch.atoms if a.atom.strip() == patch_bridge_name.strip()).number + self.offset
        anch_2_itp = anch_2_itp_obj.number
        
        stitched_mol.bonds.append(ItpBond(anch_2_itp, patch_bridge_idx, 1))

        # Remap topologies: Map N -> NZ 
        patch_n_idx = next(a.number for a in self.patch.atoms if a.atom.strip() == "N") + self.offset
        remap_dict = {
            patch_n_idx: anch_2_itp   
        }
        
        # Reindex and Sort
        reindex_topology(stitched_mol, remap_dict)
        
        # Map the PdbRecord object to the ItpAtom object
        stitched_map = self._build_pdb_to_itp_map(stitched_mol.records, stitched_mol.atoms, itp_chains=self.itp_chains)
        stitched_itp_to_pdb = {}
        for r in stitched_mol.records:
            itp_obj = stitched_map.get(id(r))
            if itp_obj:
                stitched_itp_to_pdb[itp_obj.number] = r

        # Build graph on the final, cleanly-numbered molecule
        stitched_graph = MolGraph(stitched_mol, itp_to_pdb_map=stitched_itp_to_pdb)

        # Find the PDB graph indices of the NZ and C7 atoms
        nz_idx = next(i for i, r in enumerate(stitched_mol.records) 
                    if r.res_seq == target_anchors[2].res_seq 
                    and r.name.strip() == target_anchors[2].name.strip()
                    and r.chain.strip() == target_anchors[2].chain.strip())
        c7_idx = next(i for i, r in enumerate(stitched_mol.records) 
                    if r.res_seq == target_reference.res_seq 
                    and r.name.strip() == patch_bridge_name.strip())

        # --- Translate and preserve original bonds ---
        for old_bond in self.base.bonds:
            # Look up physical PdbRecords from the old indices
            record_1 = old_itp_to_pdb.get(old_bond.a1)
            record_2 = old_itp_to_pdb.get(old_bond.a2)

            # Skip if the map missed them, or if explicitly deleted (like HZ1)
            if not record_1 or not record_2:
                continue
            if record_1 in self.base_deletions or record_2 in self.base_deletions:
                continue

            # Find new topology numbers using the stitched_map
            new_itp_1 = stitched_map.get(id(record_1))
            new_itp_2 = stitched_map.get(id(record_2))

            # If both survived, append the bond with the new numbers
            if new_itp_1 and new_itp_2:
                stitched_mol.bonds.append(ItpBond(new_itp_1.number, new_itp_2.number, old_bond.type))

        # Call checker to get the missing parameters
        angles_to_check, dihedrals_to_check = self.build_bridge_topol(stitched_mol, stitched_graph, nz_idx, c7_idx, stitched_map)
        print(f"Topology junction built: 1 bond, {len(angles_to_check)} angles, {len(dihedrals_to_check)} dihedrals, and {len(dihedrals_to_check)} pairs, added.")
        
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
        
        stitched_mol = self.update_atom_types(stitched_mol, type_map, stitched_itp_to_pdb, target_reference.chain)
        stitched_mol = self.balance_electrostatics(stitched_mol, stitched_itp_to_pdb, target_reference.chain)

        return stitched_mol

        return stitched_mol
    
    def update_atom_types(self, stitched_mol, type_map, stitched_itp_to_pdb, target_chain):
        """Updates specific atom types at the junction exclusively on the target chain."""
        for atom in stitched_mol.atoms:
            rec = stitched_itp_to_pdb.get(atom.number)
            if rec and rec.res_seq == self.res_id and rec.chain.strip() == target_chain.strip():
                if atom.atom.strip() in type_map:
                    atom.type = type_map[atom.atom.strip()]
        return stitched_mol
    
    def balance_electrostatics(self, stitched_mol, stitched_itp_to_pdb, target_chain, target_charge=0):
        """
        Balances the electrostatic charge of the patched residue to ensure a perfect 
        integer net charge. It updates the heavy anchor atoms to match the patch 
        definitions and distributes the resulting fractional deficit across specific 
        hydrogen sinks.
        
        :param stitched_mol: (Mol) The combined molecule object containing the new junction.
        :return: (Mol) The molecule with updated and balanced fractional charges.
        """
        # --- ANCHOR CHARGE INHERITANCE ---
        # The base protein anchors inherit the exact partial charges from the patch 
        # topology to ensure the new junction bond is electrostatically stable.
        anchor_map = {'NZ': 'N', 'CE': 'C10', 'CD': 'C11'}
        
        for base_atom in stitched_mol.atoms:
            rec = stitched_itp_to_pdb.get(base_atom.number)
            if rec and rec.res_seq == self.res_id and rec.chain.strip() == target_chain.strip():
                if base_atom.atom.strip() in anchor_map:
                    patch_equiv_name = anchor_map[base_atom.atom.strip()]
                    patch_equiv_atom = next(a for a in self.patch.atoms if a.atom.strip() == patch_equiv_name)
                    
                    old_charge = float(base_atom.charge)
                    base_atom.charge = float(patch_equiv_atom.charge)
                    base_atom.comment = f"; old charge: {old_charge:.4f} (anchor update)"


        # Identify which sink atoms are present in this specific residue
        possible_sinks = ['HD1' , 'HD2', 'HG1', 'HG2', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7']
        active_sinks = []
        target_residue_atoms = []

        for a in stitched_mol.atoms:
            rec = stitched_itp_to_pdb.get(a.number)
            if rec and rec.res_seq == self.res_id and rec.chain.strip() == target_chain.strip():
                target_residue_atoms.append(a)
                if a.atom.strip() in possible_sinks:
                    active_sinks.append(a)

        # Calculate the deficit
        current_sum = sum(float(a.charge) for a in target_residue_atoms)
        # target_sum = round(current_sum)
        delta_q = target_charge - current_sum     # ONLY FOR NEUTRAL LYSINES. NEED TO CHANGE THIS IF NOT.

        # Distribute only if we have sinks to take the charge
        if active_sinks and abs(delta_q) > 1e-6:
            q_per_atom = delta_q / len(active_sinks)
            for atom in active_sinks:
                old_q = float(atom.charge)
                atom.charge = old_q + q_per_atom
                atom.comment = f"; balanced (was {old_q:.4f})"

        total_stitched_charge = sum(float(a.charge) for a in target_residue_atoms)
        print(f"Total Stitched Charge for Res {self.res_id} (Chain {target_chain.strip()}): {total_stitched_charge:.6f}")
        return stitched_mol