"""
Executes the topological and geometric surgery required to attach patch molecules
to base protein targets. Handles file loading, atom deletion, index offsetting,
and electrostatic balancing.
"""

import os
from dataclasses import replace
from .geometry import MolGraph
from .topology_tools import reindex_topology
from .mol_record import Mol, ItpBond, ItpAngle, ItpDih, ItpPair
import json


class Stitcher:
    """
    Performs the topological and geometric surgery required to attach a patch molecule
    to a base protein molecule.
    """

    def __init__(self, base_mol, patch_mol, target_res_id, config, itp_chains=None):
        """
        Initializes the Stitcher and validates the configuration.

        :param base_mol: The unpatched base protein molecule.
        :type base_mol: Mol
        :param patch_mol: The patch molecule to be attached.
        :type patch_mol: Mol
        :param target_res_id: The residue sequence number of the attachment site.
        :type target_res_id: int
        :param config: The configuration dictionary for the junction parameters.
        :type config: dict
        :param itp_chains: A list of chain identifiers to restrict topology mapping, defaults to None.
        :type itp_chains: list, optional
        :raises ValueError: If the base bridge is not listed in the base anchors.
        """

        self.base = base_mol
        self.patch = patch_mol
        self.res_id = target_res_id
        self.itp_chains = itp_chains or []
        self.config = config

        if self.config["base_bridge"] not in self.config["base_anchors"]:
            raise ValueError(
                f"Config Error: 'base_bridge' ({self.config['base_bridge']}) is missing from 'base_anchors' list."
            )

        # State tracking
        self.offset = 1000000
        self.base_deletions = []
        self.patch_deletions = []
        self.junction_interactions = []

        current_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(current_dir)
        codes_path = os.path.join(project_root, "configs", "glycan_codes.json")
        with open(codes_path, "r") as f:
            self.glycan_codes = json.load(f)

    def get_bonded_hydrogens(self, mol, mol_graph, anchors):
        """
        Uses the pre-built molecular connectivity graph to identify hydrogens bonded to specific anchors.

        :param mol: The molecule object containing the coordinate records.
        :type mol: Mol
        :param mol_graph: The connectivity matrix/graph of the molecule.
        :type mol_graph: MolGraph
        :param anchors: Objects representing the anchor atoms.
        :type anchors: list
        :return: A list of PdbRecord objects representing the bonded hydrogens to be deleted.
        :rtype: list
        """

        h_to_delete = []
        for anchor in anchors:
            # Find the exact index of the anchor in the molecule's records using its serial number
            anchor_idx = next(
                i for i, r in enumerate(mol.records) if r.serial == anchor.serial
            )

            # Fetch all bonded neighbors from the networkx graph
            for n_idx in mol_graph.nx_graph.neighbors(anchor_idx):
                neighbor_record = mol.records[n_idx]
                if neighbor_record.name.strip().startswith("H"):
                    if neighbor_record not in h_to_delete:
                        h_to_delete.append(neighbor_record)
        return h_to_delete

    @staticmethod
    def _build_pdb_to_itp_map(records, atoms, itp_chains=None):
        """
        Creates a dictionary mapping the Python id() of a PdbRecord to its corresponding ItpAtom object.

        :param records: The PdbRecord objects from the coordinate file.
        :type records: list
        :param atoms: The ItpAtom objects from the topology file.
        :type atoms: list
        :param itp_chains: Chain identifiers to limit the scope of the mapping, defaults to None.
        :type itp_chains: list, optional
        :return: A dictionary bridging the physical coordinates to the topological indices.
        :rtype: dict
        """

        itp_blocks = []
        current_itp_res = None
        for a in atoms:
            if a.res_n != current_itp_res:
                itp_blocks.append({})
                current_itp_res = a.res_n
            itp_blocks[-1][a.atom.strip()] = a

        scoped_pdb_blocks = []
        current_pdb_res = None

        for r in records:
            if itp_chains and r.chain.strip() not in itp_chains:
                continue

            res_id = (r.chain.strip(), r.res_seq, r.res_name.strip())
            if res_id != current_pdb_res:
                scoped_pdb_blocks.append({})
                current_pdb_res = res_id
            scoped_pdb_blocks[-1][r.name.strip()] = r

        mapping = {}
        i, j = 0, 0

        while i < len(scoped_pdb_blocks) and j < len(itp_blocks):
            p_res = list(scoped_pdb_blocks[i].values())[0].res_name.strip()
            i_res = list(itp_blocks[j].values())[0].res.strip()

            if p_res[:2] == i_res[:2]:
                for atom_name, record_obj in scoped_pdb_blocks[i].items():
                    if atom_name in itp_blocks[j]:
                        mapping[id(record_obj)] = itp_blocks[j][atom_name]
                i += 1
                j += 1
            else:
                found_sync = False

                for offset_i in range(40):
                    if i + offset_i >= len(scoped_pdb_blocks):
                        break
                    for offset_j in range(40):
                        if j + offset_j >= len(itp_blocks):
                            break

                        sync_length = min(
                            3,
                            len(scoped_pdb_blocks) - (i + offset_i),
                            len(itp_blocks) - (j + offset_j),
                        )
                        is_match = True

                        for k in range(sync_length):
                            pk = list(scoped_pdb_blocks[i + offset_i + k].values())[
                                0
                            ].res_name.strip()[:2]
                            ik = list(itp_blocks[j + offset_j + k].values())[
                                0
                            ].res.strip()[:2]
                            if pk != ik:
                                is_match = False
                                break

                        if is_match and sync_length > 0:
                            i = i + offset_i
                            j = j + offset_j
                            found_sync = True
                            break

                    if found_sync:
                        break

                if not found_sync:
                    i += 1

        return mapping

    def build_bridge_topol(
        self,
        stitched_mol,
        stitched_graph,
        base_anchor_idx,
        patch_anchor_idx,
        stitched_map,
    ):
        """
        Uses the connectivity matrix to identify and generate missing angle and dihedral parameters.

        :param stitched_mol: The combined molecule object.
        :type stitched_mol: Mol
        :param stitched_graph: The connectivity graph of the combined molecule.
        :type stitched_graph: MolGraph
        :param base_anchor_idx: The Python list index of the base protein anchor.
        :type base_anchor_idx: int
        :param patch_anchor_idx: The Python list index of the patch anchor.
        :type patch_anchor_idx: int
        :param stitched_map: The dictionary bridging PDB memory IDs to ITP atom objects.
        :type stitched_map: dict
        :return: Two lists containing the angle tuples and dihedral tuples to be generated.
        :rtype: tuple of (list, list)
        """

        base_neighbors = list(stitched_graph.nx_graph.neighbors(base_anchor_idx))
        patch_neighbors = list(stitched_graph.nx_graph.neighbors(patch_anchor_idx))

        base_exclusive = [n for n in base_neighbors if n != patch_anchor_idx]
        patch_exclusive = [n for n in patch_neighbors if n != base_anchor_idx]

        def to_itp(pdb_idx):
            record = stitched_mol.records[pdb_idx]
            mapped_atom = stitched_map.get(id(record))

            if mapped_atom is not None:
                return mapped_atom.number

            for a in stitched_mol.atoms:
                if a.res_n == record.res_seq and a.atom.strip() == record.name.strip():
                    return a.number

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

    def _shift_topology_indices(self, top_list, offset):
        """
        Safely shifts the internal indices of GROMACS topology objects.

        :param top_list: A list of dataclass objects representing bonds, pairs, angles, or dihedrals.
        :type top_list: list
        :param offset: The integer to add to the existing indices.
        :type offset: int
        :return: A new list of shifted topology objects.
        :rtype: list
        """

        shifted = []
        for item in top_list:
            kwargs = {
                attr: getattr(item, attr) + offset
                for attr in ["a1", "a2", "a3", "a4"]
                if hasattr(item, attr)
            }
            shifted.append(replace(item, **kwargs))
        return shifted

    def stitch_molecules(self, aligned_patch_atoms, target_reference, target_anchors):
        """
        Executes the full attachment pipeline: deletes overlapping atoms, offsets topologies,
        splices lists, generates junction bonds, and balances electrostatics.

        :param aligned_patch_atoms: PdbRecord objects of the patch translated to the target site.
        :type aligned_patch_atoms: list
        :param target_reference: The primary target anchor used for chain/residue designation.
        :type target_reference: PdbRecord
        :param target_anchors: PdbRecord objects of the protein attachment site.
        :type target_anchors: list
        :raises ValueError: If the bridging anchor atom cannot be found.
        :return: A tuple containing the final stitched Mol, the new MolGraph, and the junction interaction log.
        :rtype: tuple of (Mol, MolGraph, dict)
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

        base_map = self._build_pdb_to_itp_map(
            self.base.records, self.base.atoms, itp_chains=self.itp_chains
        )

        target_res_records = [
            r
            for r in self.base.records
            if r.res_seq == self.res_id
            and r.chain.strip() == target_reference.chain.strip()
        ]
        for r in target_res_records:
            if id(r) not in base_map:
                for a in self.base.atoms:
                    if a.res_n == self.res_id and a.atmo.strip() == r.name.strip():
                        base_map[id(r)] = a
                        break

        base_bridge_record = next(
            r
            for r in self.base.records
            if r.res_seq == self.res_id
            and r.name.strip() == base_bridge_name.strip()
            and r.chain.strip() == target_reference.chain.strip()
        )

        anch_2_itp_obj = base_map.get(id(base_bridge_record))
        if anch_2_itp_obj is None:
            for a in self.base.atoms:
                if (
                    a.res_n == self.res_id
                    and a.atom.strip() == base_bridge_name.strip()
                ):
                    anch_2_itp_obj = a
                    break

        old_itp_to_pdb = {}
        for record in self.base.records:
            itp_obj = base_map.get(id(record))
            if itp_obj is not None:
                old_itp_to_pdb[itp_obj.number] = record

        self.base_deletions = [
            r
            for r in self.base.records
            if r.res_seq == self.res_id
            and r.chain.strip() == target_reference.chain.strip()
            and r.name.strip() in base_del_names
        ]
        self.patch_deletions = [
            r for r in aligned_patch_atoms if r.name.strip() in patch_del_names
        ]
        patch_overlap_anchors = [
            r for r in aligned_patch_atoms if r.name.strip() in patch_overlap_names
        ]

        print(
            f"Final kill-list: {len(self.base_deletions)} base atoms, {len(self.patch_deletions)} patch atoms."
        )

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

        max_chain_res = max(
            (
                r.res_seq
                for r in self.base.records
                if r.chain.strip() == target_reference.chain.strip()
            ),
            default=target_reference.res_seq,
        )

        for i, record in enumerate(aligned_patch_atoms):
            if (
                record not in self.patch_deletions
                and record not in patch_overlap_anchors
            ):
                original_atom = self.patch.atoms[i]

                if patch_merged_name:
                    apply_res_name = new_res_name
                    apply_res_seq = target_reference.res_seq
                else:
                    apply_res_name = record.res_name
                    apply_res_seq = max_chain_res + record.res_seq

                filter_patch_records.append(
                    replace(
                        record,
                        res_name=apply_res_name,
                        chain=target_reference.chain,
                        res_seq=apply_res_seq,
                        seg_id=dynamic_seg_id,
                    )
                )

                new_atom = replace(
                    original_atom, res_n=apply_res_seq, res=apply_res_name
                )

                if original_atom.res_n != apply_res_seq:
                    existing_comment = (
                        new_atom.comment.strip()
                        if hasattr(new_atom, "comment") and new_atom.comment
                        else ""
                    )
                    new_atom.comment = (
                        f"{existing_comment} ; old_res: {original_atom.res_n}".strip()
                    )

                filter_patch_atoms.append(new_atom)

        valid_names_str = (
            self.glycan_codes.get("res_name", {})
            .get("PDB", {})
            .get(target_res_name, target_res_name)
        )
        valid_names = valid_names_str.split()

        for r in final_records:
            if (
                r.res_seq == target_reference.res_seq
                and r.chain == target_reference.chain
            ):
                if r.res_name.strip() in valid_names:
                    r.res_name = new_res_name

        for atom in filter_patch_atoms:
            atom.number += self.offset

        patch_bonds = self._shift_topology_indices(self.patch.bonds, self.offset)
        patch_pairs = self._shift_topology_indices(self.patch.pairs, self.offset)
        patch_angles = self._shift_topology_indices(self.patch.angles, self.offset)
        patch_dihs = self._shift_topology_indices(self.patch.dihs, self.offset)

        if patch_merged_name:
            insert_idx_records = (
                next(
                    i for i, r in enumerate(final_records) if r.serial == anchor_serial
                )
                + 1
            )
        else:
            last_res_idx = max(
                i
                for i, r in enumerate(final_records)
                if r.res_seq == self.res_id
                and r.chain.strip() == target_reference.chain.strip()
            )
            insert_idx_records = last_res_idx + 1

        final_records[insert_idx_records:insert_idx_records] = filter_patch_records

        for a in final_atoms:
            rec = old_itp_to_pdb.get(a.number)
            if rec:
                new_res_n = rec.res_seq
                true_old_res = str(a.res_n)
                if "old_res:" in getattr(a, "comment", ""):
                    true_old_res = a.comment.split("old_res:")[-1].strip()

                if a.res_n != new_res_n:
                    a.res_n = new_res_n

                if str(new_res_n) == true_old_res:
                    a.comment = ""
                else:
                    a.comment = f"; old_res: {true_old_res}"

                if (
                    rec.res_seq == target_reference.res_seq
                    and rec.chain.strip() == target_reference.chain.strip()
                ):
                    if a.res.strip() == target_res_name:
                        a.res = new_res_name

        if anch_2_itp_obj is None:
            raise ValueError(
                f"Error: Anchor atom {base_bridge_name} could not be found in the topology."
            )

        if patch_merged_name:
            insert_idx_atoms = final_atoms.index(anch_2_itp_obj) + 1
        else:
            last_res_atom_idx = max(
                i
                for i, a in enumerate(final_atoms)
                if old_itp_to_pdb.get(a.number)
                and old_itp_to_pdb[a.number].res_seq == self.res_id
                and old_itp_to_pdb[a.number].chain.strip()
                == target_reference.chain.strip()
            )
            insert_idx_atoms = last_res_atom_idx + 1

        final_atoms[insert_idx_atoms:insert_idx_atoms] = filter_patch_atoms

        stitched_mol = Mol(
            self.base.name,
            self.base.moltype_section,
            final_records,
            final_atoms,
            self.base.bonds + patch_bonds,
            self.base.pairs + patch_pairs,
            self.base.angles + patch_angles,
            self.base.dihs + patch_dihs,
        )

        junction_cfg = self.config.get("itp_junction", {})
        b_type = junction_cfg.get("bond_type", 1)
        a_type = junction_cfg.get("angle_type", 5)
        d_type = junction_cfg.get("dih_type", 9)

        patch_bridge_idx = (
            next(
                a
                for a in self.patch.atoms
                if a.atom.strip() == patch_bridge_name.strip()
            ).number
            + self.offset
        )
        anch_2_itp = anch_2_itp_obj.number

        stitched_mol.bonds.append(ItpBond(anch_2_itp, patch_bridge_idx, b_type))

        remap_dict = {}
        if patch_merged_name:
            patch_n_idx = (
                next(
                    a.number
                    for a in self.patch.atoms
                    if a.atom.strip() == patch_merged_name.strip()
                )
                + self.offset
            )
            remap_dict = {patch_n_idx: anch_2_itp}

        pre_reindex_mapping = {}
        for old_num, record in old_itp_to_pdb.items():
            pre_reindex_mapping[old_num] = record

        for record, original_atom in zip(filter_patch_records, filter_patch_atoms):
            pre_reindex_mapping[original_atom.number] = record

        index_map = reindex_topology(stitched_mol, remap_dict)

        stitched_itp_to_pdb = {}
        stitched_pdb_to_itp = {}

        for old_num, record in pre_reindex_mapping.items():
            if old_num in index_map:
                new_num = index_map[old_num]
                stitched_itp_to_pdb[new_num] = record
                stitched_pdb_to_itp[id(record)] = stitched_mol.atoms[new_num - 1]

        stitched_graph = MolGraph(stitched_mol, itp_to_pdb_map=stitched_itp_to_pdb)

        if stitched_graph.is_distance_based:
            warning_msg = "Connectivity Matrix produced with distance-based method. Some glycan bonds and/or disulfide bridges may be missing."
            print(f"WARNING: {warning_msg}")
            stitched_mol.notes.append(warning_msg)

        base_idx = next(
            i
            for i, r in enumerate(stitched_mol.records)
            if r.res_seq == self.res_id
            and r.name.strip() == base_bridge_name.strip()
            and r.chain.strip() == target_reference.chain.strip()
        )

        if patch_merged_name:
            expected_patch_res_seq = target_reference.res_seq
        else:
            orig_bridge_rec = next(
                r
                for r in aligned_patch_atoms
                if r.name.strip() == patch_bridge_name.strip()
            )
            expected_patch_res_seq = max_chain_res + orig_bridge_rec.res_seq

        patch_idx = next(
            i
            for i, r in enumerate(stitched_mol.records)
            if r.res_seq == expected_patch_res_seq
            and r.name.strip() == patch_bridge_name.strip()
            and r.chain.strip() == target_reference.chain.strip()
        )

        angles_to_check, dihedrals_to_check = self.build_bridge_topol(
            stitched_mol, stitched_graph, base_idx, patch_idx, stitched_pdb_to_itp
        )
        print(
            f"Topology junction built: 1 bond, {len(angles_to_check)} angles, {len(dihedrals_to_check)} dihedrals, and {len(dihedrals_to_check)} pairs, added."
        )

        for a1, a2, a3 in angles_to_check:
            stitched_mol.angles.append(ItpAngle(a1, a2, a3, a_type))

        for a1, a2, a3, a4 in dihedrals_to_check:
            stitched_mol.dihs.append(ItpDih(a1, a2, a3, a4, d_type))
            stitched_mol.pairs.append(ItpPair(a1, a4, 1))

        stitched_mol.bonds.sort(key=lambda x: x.a1)
        stitched_mol.pairs.sort(key=lambda x: x.a1)
        stitched_mol.angles.sort(key=lambda x: x.a1)
        stitched_mol.dihs.sort(key=lambda x: x.a1)

        target_chg = self.config.get("target_charge", 0)
        stitched_mol = self.update_atom_types(
            stitched_mol,
            self.config.get("type_mapping", {}),
            stitched_itp_to_pdb,
            target_reference.chain,
        )
        stitched_mol = self.balance_electrostatics(
            stitched_mol,
            stitched_itp_to_pdb,
            target_reference.chain,
            target_charge=target_chg,
        )

        junction_log = {
            "bond": (
                index_map.get(anch_2_itp, anch_2_itp),
                index_map.get(patch_bridge_idx, patch_bridge_idx),
            ),
            "angles": angles_to_check,
            "dihs": dihedrals_to_check,
        }
        return stitched_mol, stitched_graph, junction_log

    def update_atom_types(
        self, stitched_mol, type_map, stitched_itp_to_pdb, target_chain
    ):
        """
        Updates specific atom types at the junction exclusively on the target chain.

        :param stitched_mol: The newly merged molecule object.
        :type stitched_mol: Mol
        :param type_map: A dictionary defining the old type to new type conversions.
        :type type_map: dict
        :param stitched_itp_to_pdb: Dictionary mapping ITP atoms back to PDB coordinates.
        :type stitched_itp_to_pdb: dict
        :param target_chain: The specific chain ID being modified.
        :type target_chain: str
        :return: The molecule with updated atom types.
        :rtype: Mol
        """

        self.retyped_indices = set()

        for atom in stitched_mol.atoms:
            rec = stitched_itp_to_pdb.get(atom.number)
            if (
                rec
                and rec.res_seq == self.res_id
                and rec.chain.strip() == target_chain.strip()
            ):
                if atom.atom.strip() in type_map:
                    atom.type = type_map[atom.atom.strip()]
                    self.retyped_indices.add(atom.number)

        return stitched_mol

    def balance_electrostatics(
        self, stitched_mol, stitched_itp_to_pdb, target_chain, target_charge=0
    ):
        """
        Balances the electrostatic charge of the patched residue to ensure an integer net charge.

        :param stitched_mol: The combined molecule object containing the new junction.
        :type stitched_mol: Mol
        :param stitched_itp_to_pdb: Dictionary mapping ITP atoms back to PDB coordinates.
        :type stitched_itp_to_pdb: dict
        :param target_chain: The specific chain ID being modified.
        :type target_chain: str
        :param target_charge: The desired net charge of the residue, defaults to 0.
        :type target_charge: int, optional
        :return: The molecule with updated and balanced fractional charges.
        :rtype: Mol
        """

        anchor_map = self.config.get("charge_inheritance", {})

        for base_atom in stitched_mol.atoms:
            rec = stitched_itp_to_pdb.get(base_atom.number)
            if (
                rec
                and rec.res_seq == self.res_id
                and rec.chain.strip() == target_chain.strip()
            ):
                if base_atom.atom.strip() in anchor_map:
                    patch_equiv_name = anchor_map[base_atom.atom.strip()]
                    patch_equiv_atom = next(
                        a
                        for a in self.patch.atoms
                        if a.atom.strip() == patch_equiv_name
                    )

                    old_charge = float(base_atom.charge)
                    base_atom.charge = float(patch_equiv_atom.charge)
                    base_atom.comment = (
                        f"; old charge: {old_charge:.4f} (anchor update)"
                    )

        possible_sinks = self.config.get("charge_sinks", [])
        active_sinks = []
        target_residue_atoms = []

        for a in stitched_mol.atoms:
            rec = stitched_itp_to_pdb.get(a.number)
            if (
                rec
                and rec.res_seq == self.res_id
                and rec.chain.strip() == target_chain.strip()
            ):
                target_residue_atoms.append(a)
                if a.atom.strip() in possible_sinks:
                    active_sinks.append(a)

        current_sum = sum(float(a.charge) for a in target_residue_atoms)
        delta_q = target_charge - current_sum

        if active_sinks and abs(delta_q) > 1e-6:
            q_per_atom = delta_q / len(active_sinks)
            for atom in active_sinks:
                old_q = float(atom.charge)
                atom.charge = old_q + q_per_atom
                atom.comment = f"; balanced (was {old_q:.4f})"

        total_stitched_charge = sum(float(a.charge) for a in target_residue_atoms)
        print(
            f"Total Stitched Charge for Res {self.res_id} (Chain {target_chain.strip()}): {total_stitched_charge:.6f}"
        )
        return stitched_mol

    def generate_ff_requests(self, stitched_mol, junction_log):
        """
        Translates junction atom data into corresponding atom types to request parameters.

        :param stitched_mol: The molecule containing the generated topology interactions.
        :type stitched_mol: Mol
        :param junction_log: Dictionary detailing the newly formed bonds, angles, and dihedrals.
        :type junction_log: dict
        :return: A dictionary of tuples formatted for database extraction.
        :rtype: dict
        """

        type_lookup = {atom.number: atom.type for atom in stitched_mol.atoms}
        mapped_types = set(self.config.get("type_mapping", {}).values())

        requests = {
            "bonds": set(),
            "angles": set(),
            "dihedrals_proper": set(),
            "dihedrals_improper": set(),
        }

        def is_boundary(types_tuple):
            return not all(t in mapped_types for t in types_tuple)

        if hasattr(self, "retyped_indices"):
            for b in stitched_mol.bonds:
                if b.a1 in self.retyped_indices or b.a2 in self.retyped_indices:
                    t = (type_lookup[b.a1], type_lookup[b.a2])
                    if is_boundary(t):
                        requests["bonds"].add(t)

            for a in stitched_mol.angles:
                if (
                    a.a1 in self.retyped_indices
                    or a.a2 in self.retyped_indices
                    or a.a3 in self.retyped_indices
                ):
                    t = (type_lookup[a.a1], type_lookup[a.a2], type_lookup[a.a3])
                    if is_boundary(t):
                        requests["angles"].add(t)

            for d in stitched_mol.dihs:
                if (
                    d.a1 in self.retyped_indices
                    or d.a2 in self.retyped_indices
                    or d.a3 in self.retyped_indices
                    or d.a4 in self.retyped_indices
                ):
                    t = (
                        type_lookup[d.a1],
                        type_lookup[d.a2],
                        type_lookup[d.a3],
                        type_lookup[d.a4],
                    )
                    if is_boundary(t):
                        if getattr(d, "func", 9) in [1, 9]:
                            requests["dihedrals_proper"].add(t)
                        elif getattr(d, "func", 2) in [2, 4]:
                            requests["dihedrals_improper"].add(t)

        a1, a2 = junction_log["bond"]
        requests["bonds"].add((type_lookup[a1], type_lookup[a2]))

        for a1, a2, a3 in junction_log["angles"]:
            requests["angles"].add((type_lookup[a1], type_lookup[a2], type_lookup[a3]))

        for a1, a2, a3, a4 in junction_log["dihs"]:
            requests["dihedrals_proper"].add(
                (
                    type_lookup[a1],
                    type_lookup[a2],
                    type_lookup[a3],
                    type_lookup[a4],
                )
            )

        return {k: list(v) for k, v in requests.items()}
