from mol_patcher.utilities import *
from mol_patcher.geometry import *
import networkx as nx
import json
import os

# ==========================================
# MODULE HELPER FUNCTIONS
# ==========================================

def get_atom_by_name(mol, res_seq, chain, name):
    """Fetches the PdbRecord matching name, residue, AND chain."""
    try:
        return next(r for r in mol.records 
                    if r.res_seq == res_seq 
                    and r.name.strip() == name.strip()
                    and r.chain.strip() == chain.strip())
    except StopIteration:
        print(f"Error: Could not find atom '{name}' in residue {res_seq} chain '{chain}'")
        raise

def get_downstream_atoms(nx_graph, records, pivot_idx, axis_idx):
    """Uses graph connectivity to identify the branch segment to be moved."""
    temp_nx_graph = nx_graph.copy()

    if temp_nx_graph.has_edge(pivot_idx, axis_idx):
        temp_nx_graph.remove_edge(pivot_idx, axis_idx)
    else:
        return []

    # Find the component containing the axis atom
    downstream_indices = list(nx.node_connected_component(temp_nx_graph, axis_idx))
    
    if len(downstream_indices) > 500:
        p_rec = records[pivot_idx]
        a_rec = records[axis_idx]
        raise RuntimeError(f"WARNING: Rotation at {p_rec.name}-{a_rec.name} "
                        f"attempted to move {len(downstream_indices)} atoms! "
                        "Check your graph for cycles or distance-based bonds.")

    return downstream_indices

# ==========================================
# DOMAIN-SPECIFIC SWEEPERS
# ==========================================

class RotamerSweeper:
    """Handles canonical protein sidechain rotations (Chi angles)."""
    def __init__(self, mol, mol_graph, res_seq, config, chain=" "):
        self.mol = mol
        self.graph = mol_graph
        self.res_seq = res_seq
        self.chain = chain.strip()
        self.config = config
        self.chi_definitions = self.config.get("chi_definitions", [])
        self.obj_to_idx = {id(a): i for i, a in enumerate(self.mol.records)}

    def run_sweep(self, steric_checker):
        """Tier 1 and Tier 2 Sidechain Sweeps."""
        canonical_rotamers = self.config.get("canonical_rotamers", [])
        spinner = ['|', '/', '-', '\\']

        print("Tier 1: Trying canonical staggered rotamers...")
        t1_names = self.chi_definitions[0]
        t1_pivot = self.obj_to_idx[id(get_atom_by_name(self.mol, self.res_seq, self.chain, t1_names[1]))]
        t1_axis = self.obj_to_idx[id(get_atom_by_name(self.mol, self.res_seq, self.chain, t1_names[2]))]
        t1_moving = [self.mol.records[idx] for idx in get_downstream_atoms(self.graph.nx_graph, self.mol.records, t1_pivot, t1_axis)]
        
        for i, pose in enumerate(canonical_rotamers):
            char = spinner[i % 4]
            print(f"\rTier 1 Progress... {char}", end='', flush=True)
            self.apply_pose(pose)
            
            dist, clash_info = steric_checker.check_clashes(limit_to_atoms=t1_moving)
            if dist is None:
                print(f"\rTier 1 Progress... Done. Success!{' ' * 50}")
                return True

        print(f"\rTier 1 Progress... Failed. Moving to Wiggle...{' ' * 60}")
        
        print("Tier 2: Systematic chi wiggle...")
        t2_names = self.chi_definitions[-1]
        t2_pivot = self.obj_to_idx[id(get_atom_by_name(self.mol, self.res_seq, self.chain, t2_names[1]))]
        t2_axis = self.obj_to_idx[id(get_atom_by_name(self.mol, self.res_seq, self.chain, t2_names[2]))]
        t2_moving = [self.mol.records[idx] for idx in get_downstream_atoms(self.graph.nx_graph, self.mol.records, t2_pivot, t2_axis)]

        current_pose = []
        for j in range(len(self.chi_definitions) - 1):
            names = self.chi_definitions[j]
            coords = [get_atom_by_name(self.mol, self.res_seq, self.chain, n) for n in names]
            current_pose.append(get_dihedral(*[[a.x, a.y, a.z] for a in coords]))

        for i, angle in enumerate(range(0, 360, 30)):
            char = spinner[i % 4]
            print(f"\rTier 2 Progress... {char}", end='', flush=True)
            
            self.apply_pose(current_pose + [float(angle)])
            dist, clash_info = steric_checker.check_clashes(limit_to_atoms=t2_moving)
            if dist is None:
                print(f"\rTier 2 Progress... Done. Success!{' ' * 50}")
                return True
        
        print(f"\rTier 2 Progress... Failed.{' ' * 60}")
        return False

    def apply_pose(self, chi_angles):
        for i, angle in enumerate(chi_angles):
            self.apply_chi_rotation(i, angle)

    def apply_chi_rotation(self, chi_index, target_angle):
        names = self.chi_definitions[chi_index]
        atom_records = [get_atom_by_name(self.mol, self.res_seq, self.chain, n) for n in names]
        coords = [[a.x, a.y, a.z] for a in atom_records]
        delta = target_angle - get_dihedral(*coords)
        
        pivot_idx = self.obj_to_idx[id(atom_records[1])]
        axis_idx = self.obj_to_idx[id(atom_records[2])]

        moving_indices = get_downstream_atoms(self.graph.nx_graph, self.mol.records, pivot_idx, axis_idx)
        moving_records = [self.mol.records[i] for i in moving_indices]
        rotate_dihedral(moving_records, coords[1], coords[2], delta)


class PatchSweeper:
    """Handles Tier 3 torsions for rigid patches (e.g. PFP)."""
    def __init__(self, mol, graph, res_seq, config, chain=" "):
        self.mol = mol
        self.graph = graph
        self.res_seq = res_seq
        self.chain = chain.strip()
        self.config = config
        self.obj_to_idx = {id(a): i for i, a in enumerate(self.mol.records)}

    def run_torsion(self, steric_checker):
        """Rotates only the patch molecule around the junction bond."""
        print(f"Tier 3: {self.config['target_res_name']} crowded. Twisting rigid patch junction...")
        base_bridge_name = self.config["base_bridge"]
        patch_bridge_name = self.config["patch_bridge"]

        base_atom = get_atom_by_name(self.mol, self.res_seq, self.chain, base_bridge_name)
        pivot_idx = self.obj_to_idx[id(base_atom)]
        
        patch_atom = None
        for n_idx in self.graph.nx_graph.neighbors(pivot_idx):
            if self.mol.records[n_idx].name.strip() == patch_bridge_name.strip():
                patch_atom = self.mol.records[n_idx]
                axis_idx = n_idx
                break
                
        if not patch_atom:
            print(f"Error: Could not find patch bridge '{patch_bridge_name}'.")
            return False

        moving_indices = get_downstream_atoms(self.graph.nx_graph, self.mol.records, pivot_idx, axis_idx)
        moving_records = [self.mol.records[idx] for idx in moving_indices]
        junction_records = [a for a in moving_records if a.res_seq == patch_atom.res_seq]

        spinner = ['|', '/', '-', '\\']
        for i, angle in enumerate(range(0, 360, 2)):
            char = spinner[i % 4]
            print(f"\rTier 3 Progress... {char}", end='', flush=True)

            rotate_dihedral(moving_records, (base_atom.x, base_atom.y, base_atom.z), (patch_atom.x, patch_atom.y, patch_atom.z), 2)
            dist, clash_info = steric_checker.check_clashes(limit_to_atoms=junction_records)

            if dist is None:
                print(f"\rTier 3 Progress... Done. Junction is clean!{' ' * 40}")
                return True

        print(f"\nTier 3: Could not find a completely clean junction pose.")
        return False


class GlycanSweeper:
    """Identifies and stores glycosidic linkages for optimization."""
    def __init__(self, mol, graph, library_path):
        self.mol = mol
        self.graph = graph.nx_graph
        self.library = self._load_library(library_path)
        self.linkages = []

    def _load_library(self, path):
        try:
            with open(path, 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            print(f"Warning: Glycan library not found at {path}.")
            return {}
        
    def find_linkages(self, patch_chain):
        """Scans the molecule and populates the self.linkages list."""

        iupac_map = {
            "BGLC": "Glc", "AGLC": "Glc", "GLC": "Glc", "NAG": "Glc", 
            "BMAN": "Man", "AMAN": "Man", "MAN": "Man",
            "BGAL": "Gal", "AGAL": "Gal", "GAL": "Gal",
            "AFUC": "Fuc", "BFUC": "Fuc", "FUC": "Fuc"
        }

        for node in self.graph.nodes():
            atom = self.mol.records[node]

            if atom.chain.strip() != patch_chain.strip():
                continue

            if atom.name.strip() == "C1":
                c1_idx = node
                c1_coords = (atom.x, atom.y, atom.z)
                o5_coords, c2_coords, bridge_o_idx, bridge_o_coords = None, None, None, None

                for n_idx in self.graph.neighbors(c1_idx):
                    n_atom = self.mol.records[n_idx]
                    if n_atom.name.strip() == "C2":
                        c2_idx, c2_coords = n_idx, (n_atom.x, n_atom.y, n_atom.z)
                    elif n_atom.name.strip() == "O5":
                        o5_idx, o5_coords = n_idx, (n_atom.x, n_atom.y, n_atom.z)
                    elif n_atom.name.strip().startswith("O") and n_atom.name.strip() != "O5":
                        bridge_o_idx, bridge_o_coords = n_idx, (n_atom.x, n_atom.y, n_atom.z)

                if bridge_o_idx is not None:
                    cx_idx, linkage_num = None, None
                    for nn_idx in self.graph.neighbors(bridge_o_idx):
                        nn_atom = self.mol.records[nn_idx]
                        if nn_atom.name.strip().startswith("C") and nn_idx != c1_idx:
                            cx_idx = nn_idx
                            linkage_num = nn_atom.name.strip()[1:] 
                            break 

                    if cx_idx is not None and linkage_num and linkage_num.isdigit() and c2_coords and o5_coords:
                        dihedral = get_dihedral(o5_coords, c1_coords, c2_coords, bridge_o_coords)
                        anomer = "a" if dihedral > 0 else "b"
                        
                        cx_atom = self.mol.records[cx_idx]
                        ref_num = int(linkage_num) - 1
                        ref_name = f"C{ref_num}"
                        cx1_idx, c4_idx = None, None
                        
                        for i, r in enumerate(self.mol.records):
                            if r.res_seq == cx_atom.res_seq and r.chain == cx_atom.chain:
                                if r.name.strip() == ref_name: cx1_idx = i
                                if r.name.strip() == "C4" and linkage_num == "6": c4_idx = i
                        
                        if cx1_idx is not None:
                            
                            # --- NEW: Build the Master Library Key ---
                            res1_name = atom.res_name.strip()
                            res2_name = cx_atom.res_name.strip()
                            
                            prefix1 = iupac_map.get(res1_name, "Unk")
                            prefix2 = iupac_map.get(res2_name, "Unk")
                            
                            link_type = f"{prefix1}{prefix2}_{anomer}1{linkage_num}"
                            
                            self.linkages.append({
                                "type": link_type,
                                "c1_idx": c1_idx, "ox_idx": bridge_o_idx, "cx_idx": cx_idx,
                                "o5_idx": o5_idx, "cx1_idx": cx1_idx, "c2_idx": c2_idx, "c4_idx": c4_idx
                            })

# ==========================================
# SEARCH ALGORITHMS & ORCHESTRATION
# ==========================================

class GreedyDFS:
    """Executes a Depth-First Search with backtracking across oligomer linkages."""
    def __init__(self, mol, graph, library, linkages, checker):
        self.mol = mol
        self.nx_graph = graph.nx_graph
        self.library = library
        self.linkages = linkages
        self.checker = checker

    def run_search(self):
        """Triggers the recursive search starting at the root linkage (index 0)."""
        if not self.linkages:
            return True
        return self._search(0)

    def _search(self, linkage_idx):
        """Recursive backtracking function with SNFG translation."""
        if linkage_idx == len(self.linkages):
            return True

        linkage = self.linkages[linkage_idx]
        raw_link_type = linkage["type"] # e.g., "ManGlc_b14"
        
        # --- N-GLYCAN to SNFG from Beicer ---
        nglycan_dict = {
            "GlcGlc_b14": "Nglycan_Core_01_a1",  # Core GlcNAc-GlcNAc
            "ManGlc_b14": "Nglycan_Core_02_a1",  # Core Man-GlcNAc
            "FucGlc_a16": "Nglycan_Core_03_a1",  # Core Fuc-GlcNAc
            "ManMan_a13": "Nglycan_Core_04_a1",  # Man-Man a1-3 (Branch)
            "ManMan_a16": "Nglycan_Core_05_a1",  # Man-Man a1-6 (Branch)
            "ManMan_a12": "Nglycan_Core_06_a1",  # Man-Man a1-2 (Branch)
            "GalGlc_b14": "Nglycan_Core_09_a1",  # Gal-GlcNAc
            "GlcMan_b12": "Nglycan_Core_10_a1",  # GlcNAc-Man b1-2
            "GlcMan_b14": "Nglycan_Core_12_a1",  # GlcNAc-Man b1-4
            "GlcMan_b16": "Nglycan_Core_14_a1",  # GlcNAc-Man b1-6
            "FucGlc_a13": "Nglycan_Core_15_a1"   # Fuc-GlcNAc a1-3
        }

        # Try specific N-glycan data
        if raw_link_type in nglycan_dict:
            target_key = nglycan_dict[raw_link_type]
            basins = self.library.get(target_key, [])
        else:
            basins = []

        # Fall back to generic homodimers (e.g., standard ManMan_a13)
        if not basins:
            basins = self.library.get(raw_link_type, [])
            
        # Final physical proxy (Forces generic GlcGlc geometry)
        if not basins:
            bond_type = raw_link_type.split("_")[-1] if "_" in raw_link_type else raw_link_type
            fallback_key = f"GlcGlc_{bond_type}"
            basins = self.library.get(fallback_key, [])
            
            if not basins:
                print(f"Warning: {raw_link_type} completely missing. Leaving frozen.")
                return self._search(linkage_idx + 1)

        c1_idx, ox_idx, cx_idx = linkage["c1_idx"], linkage["ox_idx"], linkage["cx_idx"]
        cx1_idx, c2_idx, c4_idx = linkage["cx1_idx"], linkage["c2_idx"], linkage["c4_idx"]

        # --- REVERSED AXIS FIX ---
        downstream_phi = get_downstream_atoms(self.nx_graph, self.mol.records, pivot_idx=ox_idx, axis_idx=c1_idx)
        downstream_psi = get_downstream_atoms(self.nx_graph, self.mol.records, pivot_idx=cx_idx, axis_idx=ox_idx)
        
        # Combine your raw index lists first
        all_moving_indices = downstream_phi + downstream_psi

        # Map the indices (integers) to their current coordinates
        restore_point = {idx: (self.mol.records[idx].x, self.mol.records[idx].y, self.mol.records[idx].z) 
                        for idx in all_moving_indices}
        
        # Create object list for rotations/sterics afterwards
        moving_records = [self.mol.records[idx] for idx in all_moving_indices]
        
        for basin in basins:
            c1_coords = [self.mol.records[c1_idx].x, self.mol.records[c1_idx].y, self.mol.records[c1_idx].z]
            c2_coords = [self.mol.records[c2_idx].x, self.mol.records[c2_idx].y, self.mol.records[c2_idx].z]
            ox_coords = [self.mol.records[ox_idx].x, self.mol.records[ox_idx].y, self.mol.records[ox_idx].z]
            cx_coords = [self.mol.records[cx_idx].x, self.mol.records[cx_idx].y, self.mol.records[cx_idx].z]
            cx1_coords = [self.mol.records[cx1_idx].x, self.mol.records[cx1_idx].y, self.mol.records[cx1_idx].z]

            delta_phi = basin["phi"] - get_dihedral(c2_coords, c1_coords, ox_coords, cx_coords)
            delta_psi = basin["psi"] - get_dihedral(c1_coords, ox_coords, cx_coords, cx1_coords)

            rotate_dihedral([self.mol.records[i] for i in downstream_phi], c1_coords, ox_coords, delta_phi)
            rotate_dihedral([self.mol.records[i] for i in downstream_psi], ox_coords, cx_coords, delta_psi)

            moving_records = [self.mol.records[i] for i in downstream_phi + downstream_psi]

            if basin.get("omega") is not None and c4_idx is not None:
                c4_coords = [self.mol.records[c4_idx].x, self.mol.records[c4_idx].y, self.mol.records[c4_idx].z]
                delta_omega = basin["omega"] - get_dihedral(ox_coords, cx_coords, cx1_coords, c4_coords)
                
                downstream_omega = get_downstream_atoms(self.nx_graph, self.mol.records, pivot_idx=cx1_idx, axis_idx=cx_idx)
                rotate_dihedral([self.mol.records[i] for i in downstream_omega], cx1_coords, cx_coords, delta_omega)
                moving_records.extend([self.mol.records[i] for i in downstream_omega])

            child_res_seq = self.mol.records[c1_idx].res_seq
            local_moving = [a for a in moving_records if a.res_seq == child_res_seq]

            dist, msg = self.checker.check_clashes(limit_to_atoms=local_moving)
            
            if dist is None:
                if self._search(linkage_idx + 1):
                    return True
            else:
                # --- ADD THIS BACK IN ---
                print(f"      [Debug] Clash on {raw_link_type} (Basin {basin['phi']}/{basin['psi']}): {msg}")
            
            for i, (old_x, old_y, old_z) in restore_point.items():
                self.mol.records[i].x = old_x
                self.mol.records[i].y = old_y
                self.mol.records[i].z = old_z
                    
        return False


class SweepConductor:
    """Orchestrates the various sweeps for steric resolution."""
    def __init__(self, mol, graph, target_resid, config, project_root, chain_id):
        self.mol = mol
        self.graph = graph
        self.resid = target_resid
        self.config = config
        self.chain = chain_id.strip()
        
        self.prot_sweeper = RotamerSweeper(mol, graph, target_resid, config, chain=self.chain)
        self.is_oligomer = "patch_merged_atom" not in config
        
        if self.is_oligomer:
            lib_path = os.path.join(project_root, 'configs', 'glycosidic_angles.json')
            self.glycan_sweeper = GlycanSweeper(mol, graph, lib_path)
            self.glycan_sweeper.find_linkages(patch_chain=self.chain)
        else:
            self.patch_sweeper = PatchSweeper(mol, graph, target_resid, config, chain=self.chain)

    def optimize(self, moving_atoms, static_atoms):
        checker = StericChecker(self.graph, moving_atoms, static_atoms, tolerance=0.75)
        
        if self.is_oligomer:
            print("\nStarting Global Root-to-Leaf Oligomer Sweep...")
            canonical_rotamers = self.config.get("canonical_rotamers", [])
            sidechain_atoms = [a for a in moving_atoms if a.res_seq == self.resid and a.chain.strip() == self.chain]

            dfs_engine = GreedyDFS(self.mol, self.graph, self.glycan_sweeper.library, self.glycan_sweeper.linkages, checker)

            for i, pose in enumerate(canonical_rotamers):
                print(f"\n--- Testing ASN Pose {i+1}: {pose} ---")
                self.prot_sweeper.apply_pose(pose)
                dist, msg = checker.check_clashes(limit_to_atoms=sidechain_atoms)
                
                if dist is None:
                    print(f"ASN Pose Safe! Starting Glycan DFS...")
                    if dfs_engine.run_search():
                        print("\nSUCCESS: Found clash-free global pose.")
                        return True
                    else:
                        print(f"Glycan DFS failed for ASN pose {pose}. Backtracking...")
                else:
                    print(f"ASN Pose Rejected: {msg}")
            
            print("\nFAILURE: Dead end reached. All poses exhausted.")
            return False
            
        else:
            print("\nStarting Rigid Patch Optimization...")
            if not self.prot_sweeper.run_sweep(checker):
                return self.patch_sweeper.run_torsion(checker)
            return True