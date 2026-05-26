from mol_patcher.utilities import get_dihedral
from mol_patcher.geometry import StericChecker, rotate_dihedral, random_rot
import networkx as nx
import numpy as np


def get_atom_by_name(mol, res_seq, chain, name):
    """Fetches the PdbRecord matching name, residue, AND chain."""
    try:
        return next(
            r
            for r in mol.records
            if r.res_seq == res_seq
            and r.name.strip() == name.strip()
            and r.chain.strip() == chain.strip()
        )
    except StopIteration:
        print(
            f"Error: Could not find atom '{name}' in residue {res_seq} chain '{chain}'"
        )
        raise


def get_downstream_atoms(nx_graph, records, pivot_idx, axis_idx):
    """Uses graph connectivity to identify the branch segment to be moved."""
    temp_nx_graph = nx_graph.copy()

    if temp_nx_graph.has_edge(pivot_idx, axis_idx):
        temp_nx_graph.remove_edge(pivot_idx, axis_idx)
    else:
        return []

    return list(nx.node_connected_component(temp_nx_graph, axis_idx))


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

    def apply_pose(self, chi_angles):
        for i, angle in enumerate(chi_angles):
            self.apply_chi_rotation(i, angle)

    def apply_chi_rotation(self, chi_index, target_angle):
        names = self.chi_definitions[chi_index]
        atom_records = [
            get_atom_by_name(self.mol, self.res_seq, self.chain, n) for n in names
        ]
        coords = [[a.x, a.y, a.z] for a in atom_records]
        delta = target_angle - get_dihedral(*coords)

        pivot_idx = self.obj_to_idx[id(atom_records[1])]
        axis_idx = self.obj_to_idx[id(atom_records[2])]

        moving_indices = get_downstream_atoms(
            self.graph.nx_graph, self.mol.records, pivot_idx, axis_idx
        )
        moving_records = [self.mol.records[i] for i in moving_indices]
        rotate_dihedral(moving_records, coords[1], coords[2], delta)

    def get_best_pose(self, steric_checker, moving_atoms):
        """Tests all canonical sidechain poses, scores them, and returns the best state."""
        canonical_rotamers = self.config.get("canonical_rotamers", [])

        best_penalty = float("inf")
        best_pose = None
        best_coords = None

        # Save initial state of the moving cluster using indices
        base_coords = {self.obj_to_idx[id(a)]: (a.x, a.y, a.z) for a in moving_atoms}

        for i, pose in enumerate(canonical_rotamers):
            self.apply_pose(pose)

            # Unpack the tuple correctly
            penalty, count = steric_checker.score_pose(limit_to_atoms=moving_atoms)

            if penalty == 0:
                print(f"      Perfect fit found at Pose {pose}")
                return (
                    penalty,
                    pose,
                    {self.obj_to_idx[id(a)]: (a.x, a.y, a.z) for a in moving_atoms},
                )

            if penalty < best_penalty:
                best_penalty = penalty
                best_pose = pose
                best_coords = {
                    self.obj_to_idx[id(a)]: (a.x, a.y, a.z) for a in moving_atoms
                }

            # Restore coordinates for the next test
            for idx, (x, y, z) in base_coords.items():
                self.mol.records[idx].x = x
                self.mol.records[idx].y = y
                self.mol.records[idx].z = z

        return best_penalty, best_pose, best_coords


class SCRSweeper:
    """Docking engine for non-covalently bound SCRs via Monte Carlo Sampling"""

    def __init__(self, mol, graph, config, checker):
        self.mol = mol
        self.graph = graph
        self.config = config
        self.checker = checker

    def dock_near(self, glycan_atoms, scr_atoms, max_dist_nm=1.5):

        if not glycan_atoms or not scr_atoms:
            print("Error: Atom lists are empty. Cannot calculate centroid.")
            return False

        checker = StericChecker(self.graph, scr_atoms.records, self.mol.records)
        # Glycan Centroid
        glycan_coords = np.array([[a.x, a.y, a.z] for a in glycan_atoms])
        centroid = np.mean(glycan_coords, axis=0)

        base_scr_coords = np.array([[a.x, a.y, a.z] for a in scr_atoms.records])

        # Random Sampling (Monte Carlo)
        best_penalty = float("inf")
        best_pose_coords = None

        print(
            f"Docking SCR near Glycan Centroid [{centroid[0]:.2f}, {centroid[1]:.2f}, {centroid[2]:.2f}]..."
        )

        for _ in range(500):  # 500 attempts
            temp_coords = random_rot(
                scr_atoms.records, centroid, base_scr_coords, max_dist_nm
            )
            # Score
            penalty, _ = checker.score_pose(limit_to_atoms=scr_atoms.records)

            if penalty < best_penalty:
                best_penalty = penalty
                best_pose_coords = temp_coords

        # Apply best
        if best_pose_coords is not None:
            for i, atom in enumerate(scr_atoms.records):
                atom.x, atom.y, atom.z = best_pose_coords[i]

            print(f"SCR Docking complete. Best penalty: {best_penalty:.2f}")
            return True
        else:
            print("Error: Could not calculate valid poses")
            return False


class SweepConductor:
    """Orchestrates the sweeps for steric resolution."""

    def __init__(self, mol, graph, target_resid, config, chain_id):
        self.mol = mol
        self.graph = graph
        self.resid = target_resid
        self.config = config
        self.chain = chain_id.strip()

        self.prot_sweeper = RotamerSweeper(
            mol, graph, target_resid, config, chain=self.chain
        )

    def optimize(self, moving_atoms, static_atoms):
        checker = StericChecker(self.graph, moving_atoms, static_atoms)

        print("\nStarting Sweep...")
        penalty, pose, coords = self.prot_sweeper.get_best_pose(checker, moving_atoms)

        if coords is not None:
            print("\nSUCCESS: Applied best global configuration.")
            print(f"  -> Chi Pose  : {pose}")
            print(f"  -> Collision Penalty: {penalty:.2f}")

            # Apply the winning coordinates to the molecule
            for idx, (x, y, z) in coords.items():
                self.mol.records[idx].x = x
                self.mol.records[idx].y = y
                self.mol.records[idx].z = z
            return True
        else:
            print("\nFAILURE: Could not calculate any poses.")
            return False
