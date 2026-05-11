import os
import numpy as np
from dataclasses import replace
from .geometry import MolGraph
from .utilities import get_distance
from .topology_tools import reindex_topology
from .mol_record import Mol, ItpBond, ItpAngle, ItpDih, ItpPair, ItpAtom
from .pdb_io import PdbParser
import json

class Patchloader:
    """
    Locates and loads the standard patch ligand PDB and ITP files into Mol objects for processing.
    """

    def __init__(self, config_file):
        self.f = os.path.abspath(__file__)
        self.mol_patcher_dir = os.path.dirname(self.f)
        self.root = os.path.dirname(self.mol_patcher_dir)
        self.pdb_dir = os.path.join(self.root, "pdbs")
        self.top_dir = os.path.join(self.root, "topology_files")
        self.configs = os.path.join(self.root, 'configs', config_file)
        

        with open(self.configs, 'r') as file:
            self.setup = json.load(file)
    
    def get_pfp_pdb(self):
        patch_pdb_name = self.setup['patch_files']['pdb']
        patch_pdb_path = os.path.join(self.pdb_dir, patch_pdb_name)
        parser = PdbParser()
        _, patch_atoms, _ = parser.read_file(patch_pdb_path)
        print(f'patch pdb: {patch_pdb_path}')
        return patch_atoms

    
    def get_pfp_itp(self, mol_obj):
        # Loads the ITP topology into an existing Mol object.
        patch_itp_name = self.setup['patch_files']['itp']
        patch_itp_path = os.path.join(self.top_dir, patch_itp_name)
        mol_obj.load_itp(patch_itp_path)
        print(f'Patch itp: {patch_itp_path}')


class Stitcher:
    """
    Performs the topological and geometric surgery required to attach a patch molecule 
    to a base protein molecule.
    """
    def __init__(self, base_mol, patch_mol, target_res_id, config, itp_chains=None):
        self.base = base_mol
        self.patch = patch_mol
        self.res_id = target_res_id
        self.itp_chains = itp_chains or []
        self.config = config

        if self.config["base_bridge"] not in self.config["base_anchors"]:
            raise ValueError(f"Config Error: 'base_bridge' ({self.config['base_bridge']}) is missing from 'base_anchors' list.")
        
        # State tracking
        self.offset = 100000
        self.base_deletions = []
        self.patch_deletions = []
        self.junction_interactions = []

        current_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(current_dir)
        codes_path = os.path.join(project_root, 'configs', 'glycan_codes.json')
        with open(codes_path, 'r') as f:
            self.glycan_codes = json.load(f)

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

    def build_bridge_topol(self, stitched_mol, stitched_graph, base_anchor_idx, patch_anchor_idx, stitched_map):
        """
        Uses the connectivity matrix to identify and generate missing angle and dihedral 
        topology parameters across the newly formed junction bond.
        
        :param stitched_mol: (Mol) The combined molecule object.
        :param stitched_graph: (MolGraph) The connectivity graph of the combined molecule.
        :param base_anchor_idx: (int) The Python list index of the base protein anchor.
        :param patch_anchor_idx: (int) The Python list index of the patch anchor.
        :param stitched_map: (dict) The map relating the pdb indices to the itp atom indices
        :return: (tuple) Two lists containing the angle tuples and dihedral tuples to be appended.
        """
        
        base_neighbors = list(stitched_graph.nx_graph.neighbors(base_anchor_idx))
        patch_neighbors = list(stitched_graph.nx_graph.neighbors(patch_anchor_idx))
        
        base_exclusive = [n for n in base_neighbors if n != patch_anchor_idx]
        patch_exclusive = [n for n in patch_neighbors if n != base_anchor_idx]

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
        base_itp = to_itp(base_anchor_idx)
        patch_itp = to_itp(patch_anchor_idx)

        angles_to_check = []
        for n in base_exclusive:
            angles_to_check.append((to_itp(n), base_itp, patch_itp))
        for c in patch_exclusive:
            angles_to_check.append((base_itp, patch_itp, to_itp(c)))

        dihedrals_to_check = []
        for n in base_exclusive:
            for c in patch_exclusive:
                dihedrals_to_check.append((to_itp(n), base_itp, patch_itp, to_itp(c)))

        return angles_to_check, dihedrals_to_check
    

    def stitch_molecules(self, aligned_patch_atoms, target_reference, target_anchors):
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
        
        # VARIABLES (JSON)
        base_del_names = self.config["base_deletions"]
        patch_overlap_names = self.config["patch_anchors"]
        patch_del_names = self.config.get("patch_deletions", []) 
        
        base_bridge_name = self.config["base_bridge"]
        patch_bridge_name = self.config["patch_bridge"]
        patch_merged_name = self.config.get("patch_merged_atom", None)
        
        target_res_name = self.config["target_res_name"]
        new_res_name = self.config["new_res_name"]

        anchor_serial = target_anchors[2].serial
        anchor_res_seq = target_anchors[2].res_seq
        anchor_name = target_anchors[2].name.strip()

        # MAPPING AND DELETION LISTS
        base_map = self._build_pdb_to_itp_map(self.base.records, self.base.atoms, itp_chains=self.itp_chains)

        base_bridge_record = next(r for r in self.base.records 
                                if r.res_seq == self.res_id 
                                and r.name.strip() == base_bridge_name.strip()
                                and r.chain.strip() == target_reference.chain.strip())

        anch_2_itp_obj = base_map.get(id(base_bridge_record))
        if anch_2_itp_obj is None:
            for a in self.base.atoms:
                if a.res_n == self.res_id and a.atom.strip() == base_bridge_name.strip():
                    anch_2_itp_obj = a
                    break
                
        old_itp_to_pdb = {}
        for record in self.base.records:
            itp_obj = base_map.get(id(record))
            if itp_obj is not None:
                old_itp_to_pdb[itp_obj.number] = record

        self.base_deletions = [r for r in self.base.records 
                            if r.res_seq == self.res_id 
                            and r.chain.strip() == target_reference.chain.strip() 
                            and r.name.strip() in base_del_names]
        self.patch_deletions = [r for r in aligned_patch_atoms if r.name.strip() in patch_del_names]
        patch_overlap_anchors = [r for r in aligned_patch_atoms if r.name.strip() in patch_overlap_names]

        print(f"Final kill-list: {len(self.base_deletions)} base atoms, {len(self.patch_deletions)} patch atoms.")

        # FILTER AND RENAME
        final_records = [r for r in self.base.records if r not in self.base_deletions]
        
        final_atoms = []
        for a in self.base.atoms:
            rec = old_itp_to_pdb.get(a.number)
            if rec in self.base_deletions:
                continue
            final_atoms.append(a)
        
        filter_patch_records, filter_patch_atoms = [], []
        dynamic_seg_id = target_anchors[0].seg_id

        # Append patch atoms with their new designation
        # Find the highest residue number currently in the target chain if for new residues from the patch
        max_chain_res = max((r.res_seq for r in self.base.records if r.chain.strip() == target_reference.chain.strip()), default=target_reference.res_seq)

        # Append patch atoms with their new designation
        for i, record in enumerate(aligned_patch_atoms):
            if record not in self.patch_deletions and record not in patch_overlap_anchors:
                original_atom = self.patch.atoms[i]
                
                # If merging into one residue (PFP), squash the name and number to the anchor site.
                # If appending a glycan, keep the sugar name and add it to the END of the chain sequence.
                if patch_merged_name:
                    apply_res_name = new_res_name
                    apply_res_seq = target_reference.res_seq
                else:
                    apply_res_name = record.res_name
                    apply_res_seq = max_chain_res + record.res_seq 

                filter_patch_records.append(replace(record, res_name=apply_res_name, 
                                        chain=target_reference.chain, res_seq=apply_res_seq, 
                                        seg_id=dynamic_seg_id))
                
                new_atom = replace(original_atom, res_n=apply_res_seq, res=apply_res_name)

                if original_atom.res_n != apply_res_seq:
                    existing_comment = new_atom.comment.strip() if hasattr(new_atom, 'comment') and new_atom.comment else ""
                    new_atom.comment = f"{existing_comment} ; old_res: {original_atom.res_n}".strip()
                    
                filter_patch_atoms.append(new_atom)

        valid_names_str = self.glycan_codes.get("res_name", {}).get("PDB", {}).get(target_res_name, target_res_name)
        valid_names = valid_names_str.split()

        # Rename the preserved base protein atoms at the target site
        for r in final_records:
            if r.res_seq == target_reference.res_seq and r.chain == target_reference.chain:
                if r.res_name.strip() in valid_names:
                    r.res_name = new_res_name

        for atom in filter_patch_atoms:
            atom.number += self.offset

        patch_bonds = [replace(b, a1=b.a1+self.offset, a2=b.a2+self.offset) for b in self.patch.bonds]
        patch_pairs = [replace(p, a1=p.a1+self.offset, a2=p.a2+self.offset) for p in self.patch.pairs]
        patch_angles = [replace(ang, a1=ang.a1+self.offset, a2=ang.a2+self.offset, a3=ang.a3+self.offset) for ang in self.patch.angles]
        patch_dihs = [replace(d, a1=d.a1+self.offset, a2=d.a2+self.offset, a3=d.a3+self.offset, a4=d.a4+self.offset) for d in self.patch.dihs]

        # Splicing and ITP sync
        if patch_merged_name:
            insert_idx_records = next(i for i, r in enumerate(final_records) if r.serial == anchor_serial) + 1
        else:
            # CHANGED: Finds the max index of the TARGET RESIDUE instead of the whole chain
            last_res_idx = max(i for i, r in enumerate(final_records) if r.res_seq == self.res_id and r.chain.strip() == target_reference.chain.strip())
            insert_idx_records = last_res_idx + 1
            
        final_records[insert_idx_records:insert_idx_records] = filter_patch_records
        
        # Sync base ITP residue numbers to the PDB and track old numbers
        for a in final_atoms:
            rec = old_itp_to_pdb.get(a.number)
            if rec:
                new_res_n = rec.res_seq
                true_old_res = str(a.res_n)
                if 'old_res:' in getattr(a, 'comment', ''):
                    true_old_res = a.comment.split('old_res:')[-1].strip()
                    
                if a.res_n != new_res_n:
                    a.res_n = new_res_n
                    
                if str(new_res_n) == true_old_res:
                    a.comment = ""
                else:
                    a.comment = f"; old_res: {true_old_res}"
            
                if rec.res_seq == target_reference.res_seq and rec.chain.strip() == target_reference.chain.strip():
                    if a.res.strip() == target_res_name:
                        a.res = new_res_name
        
        if anch_2_itp_obj is None:
            raise ValueError(f"Error: Anchor atom {base_bridge_name} could not be found in the topology.")

        if patch_merged_name:
            insert_idx_atoms = final_atoms.index(anch_2_itp_obj) + 1
        else:
            last_res_atom_idx = max(
                i for i, a in enumerate(final_atoms) 
                if old_itp_to_pdb.get(a.number) 
                and old_itp_to_pdb[a.number].res_seq == self.res_id
                and old_itp_to_pdb[a.number].chain.strip() == target_reference.chain.strip()
            )
            insert_idx_atoms = last_res_atom_idx + 1

        final_atoms[insert_idx_atoms:insert_idx_atoms] = filter_patch_atoms

        stitched_mol = Mol(self.base.name, self.base.moltype_section, final_records, final_atoms, 
                        self.base.bonds + patch_bonds, self.base.pairs + patch_pairs, 
                        self.base.angles + patch_angles, self.base.dihs + patch_dihs)
        
        # JUNCTION INTERACTIONS AND TOPOLOGY REINDEXING
        junction_cfg = self.config.get("itp_junction", {})
        b_type = junction_cfg.get("bond_type", 1)
        a_type = junction_cfg.get("angle_type", 5)
        d_type = junction_cfg.get("dih_type", 9)

        patch_bridge_idx = next(a for a in self.patch.atoms if a.atom.strip() == patch_bridge_name.strip()).number + self.offset
        anch_2_itp = anch_2_itp_obj.number

        stitched_mol.bonds.append(ItpBond(anch_2_itp, patch_bridge_idx, b_type))

        # Initialize the remap dictionary
        remap_dict = {}
        if patch_merged_name:
            patch_n_idx = next(a.number for a in self.patch.atoms if a.atom.strip() == patch_merged_name.strip()) + self.offset
            remap_dict = {patch_n_idx: anch_2_itp}
            
        # Capture pre-reindex mapping
        pre_reindex_mapping = {}
        for old_num, record in old_itp_to_pdb.items():
            pre_reindex_mapping[old_num] = record
            
        for record, original_atom in zip(filter_patch_records, filter_patch_atoms):
            # original_atom.number still has the 100000 offset here
            pre_reindex_mapping[original_atom.number] = record

        # REINDEX THE TOPOLOGY (Mods atom.number in place)
        index_map = reindex_topology(stitched_mol, remap_dict)

        # BUILD LOOKUP DICTIONARIES
        stitched_itp_to_pdb = {}
        stitched_pdb_to_itp = {}
        
        for old_num, record in pre_reindex_mapping.items():
            if old_num in index_map:
                new_num = index_map[old_num]
                stitched_itp_to_pdb[new_num] = record
                stitched_pdb_to_itp[id(record)] = stitched_mol.atoms[new_num - 1]

        # BUILD THE GRAPH
        stitched_graph = MolGraph(stitched_mol, itp_to_pdb_map=stitched_itp_to_pdb)

        if stitched_graph.is_distance_based:
            warning_msg = "Connectivity Matrix produced with distance-based method. Some glycan bonds and/or disulfide bridges may be missing."
            print(f"WARNING: {warning_msg}")
            stitched_mol.notes.append(warning_msg)

        base_idx = next(i for i, r in enumerate(stitched_mol.records)
                        if r.res_seq == self.res_id
                        and r.name.strip() == base_bridge_name.strip()
                        and r.chain.strip() == target_reference.chain.strip())
        
        if patch_merged_name:
            expected_patch_res_seq = target_reference.res_seq
        else:
            orig_bridge_rec = next(r for r in aligned_patch_atoms if r.name.strip() == patch_bridge_name.strip())
            expected_patch_res_seq = max_chain_res + orig_bridge_rec.res_seq

        patch_idx = next(i for i, r in enumerate(stitched_mol.records)
                        if r.res_seq == expected_patch_res_seq
                        and r.name.strip() == patch_bridge_name.strip()
                        and r.chain.strip() == target_reference.chain.strip())

        # GENERATE BRIDGE TOPOLOGY 
        # (Pass stitched_pdb_to_itp to ensure the dihedral builder grabs the correct atoms)
        angles_to_check, dihedrals_to_check = self.build_bridge_topol(stitched_mol, stitched_graph, base_idx, patch_idx, stitched_pdb_to_itp)
        print(f"Topology junction built: 1 bond, {len(angles_to_check)} angles, {len(dihedrals_to_check)} dihedrals, and {len(dihedrals_to_check)} pairs, added.")
        
        for a1, a2, a3 in angles_to_check:
            stitched_mol.angles.append(ItpAngle(a1, a2, a3, a_type))

        for a1, a2, a3, a4 in dihedrals_to_check:
            stitched_mol.dihs.append(ItpDih(a1, a2, a3, a4, d_type))
            stitched_mol.pairs.append(ItpPair(a1, a4, 1))
            
        stitched_mol.bonds.sort(key=lambda x: x.a1)
        stitched_mol.pairs.sort(key=lambda x: x.a1)
        stitched_mol.angles.sort(key=lambda x: x.a1)
        stitched_mol.dihs.sort(key=lambda x: x.a1)

        # ELECTROSTATICS & TYPING
        target_chg = self.config.get("target_charge", 0) 
        stitched_mol = self.update_atom_types(stitched_mol, self.config.get("type_mapping", {}), stitched_itp_to_pdb, target_reference.chain)
        stitched_mol = self.balance_electrostatics(stitched_mol, stitched_itp_to_pdb, target_reference.chain, target_charge=target_chg)

        junction_log = {
            "bond":(index_map.get(anch_2_itp, anch_2_itp), index_map.get(patch_bridge_idx, patch_bridge_idx)),
            "angles": angles_to_check,
            "dihs": dihedrals_to_check
        }
        return stitched_mol, stitched_graph, junction_log
    
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
        # ANCHOR CHARGE INHERITANCE
        # The base protein anchors inherit the exact partial charges from the patch 
        # topology to ensure the new junction bond is electrostatically stable.
        anchor_map = self.config.get("charge_inheritance", {})
        
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
        possible_sinks = self.config.get("charge_sinks", [])
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
