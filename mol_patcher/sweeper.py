"""
Orchestrates the conformational optimization of flexible sidechains and payloads.
Uses graph traversal to isolate moving branches and sweeps canonical rotamer states
to resolve steric clashes.
"""

from mol_patcher.utilities import get_dihedral
from mol_patcher.geometry import StericChecker, rotate_dihedral
import networkx as nx

STANDARD_AMINO_ACIDS = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLU",
    "GLN",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "HSD",
    "HSE",
    "HSP",
    "GLH",
    "ASH",
    "CYX",
}


def get_atom_by_name(mol, res_seq, chain, name):
    """
    Fetches a specific atom record based on its sequence, chain, and name.

    :param mol: The molecule object containing the records.
    :type mol: Mol
    :param res_seq: The residue sequence number.
    :type res_seq: int
    :param chain: The single-character chain identifier.
    :type chain: str
    :param name: The atom name (e.g., 'CA', 'CB').
    :type name: str
    :raises StopIteration: If the specified atom cannot be found in the molecule.
    :return: The matching atom record.
    :rtype: PdbRecord
    """

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
    """
    Uses graph connectivity to isolate the downstream segment of a molecule that
    will move when a specific bond is rotated.

    :param nx_graph: The molecular connectivity graph.
    :type nx_graph: networkx.Graph
    :param records: The full list of atom records.
    :type records: list
    :param pivot_idx: The list index of the stationary atom.
    :type pivot_idx: int
    :param axis_idx: The list index of the moving base atom.
    :type axis_idx: int
    :return: A list of indices representing the isolated downstream fragment.
    :rtype: list
    """

    temp_nx_graph = nx_graph.copy()

    if temp_nx_graph.has_edge(pivot_idx, axis_idx):
        temp_nx_graph.remove_edge(pivot_idx, axis_idx)
    else:
        return []

    return list(nx.node_connected_component(temp_nx_graph, axis_idx))


class RotamerSweeper:
    """Handles canonical protein sidechain rotations (Chi angles)."""

    def __init__(self, mol, mol_graph, res_seq, config, chain=" "):
        """
        Initializes the RotamerSweeper.

        :param mol: The combined molecule object.
        :type mol: Mol
        :param mol_graph: The connectivity graph of the molecule.
        :type mol_graph: MolGraph
        :param res_seq: The target residue sequence number.
        :type res_seq: int
        :param config: The configuration dictionary containing chi definitions.
        :type config: dict
        :param chain: The target chain identifier, defaults to " ".
        :type chain: str, optional
        """

        self.mol = mol
        self.graph = mol_graph
        self.res_seq = res_seq
        self.chain = chain.strip()
        self.config = config
        self.chi_definitions = self.config.get("chi_definitions", [])
        self.obj_to_idx = {id(a): i for i, a in enumerate(self.mol.records)}

    def apply_pose(self, chi_angles):
        """
        Iterates through a list of target chi angles and applies them sequentially.

        :param chi_angles: A list of target angles in degrees for each defined chi junction.
        :type chi_angles: list
        """

        for i, angle in enumerate(chi_angles):
            self.apply_chi_rotation(i, angle)

    def apply_chi_rotation(self, chi_index, target_angle):
        """
        Calculates the delta required to reach a target angle and applies the rotation.

        :param chi_index: The index of the chi definition in the config array.
        :type chi_index: int
        :param target_angle: The desired dihedral angle in degrees.
        :type target_angle: float
        """

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
        """
        Tests all canonical sidechain poses, evaluates them using the StericChecker,
        and returns the best non-clashing state.

        :param steric_checker: An initialized StericChecker object.
        :type steric_checker: StericChecker
        :param moving_atoms: The full list of atoms considered part of the moving branch.
        :type moving_atoms: list
        :return: A tuple containing the best penalty score, the optimal pose array,
                 and a dictionary of the winning coordinates.
        :rtype: tuple of (float, list, dict)
        """

        canonical_rotamers = self.config.get("canonical_rotamers", [])

        best_penalty = float("inf")
        best_pose = None
        best_coords = None

        base_coords = {self.obj_to_idx[id(a)]: (a.x, a.y, a.z) for a in moving_atoms}

        for i, pose in enumerate(canonical_rotamers):
            self.apply_pose(pose)

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

            for idx, (x, y, z) in base_coords.items():
                self.mol.records[idx].x = x
                self.mol.records[idx].y = y
                self.mol.records[idx].z = z

        return best_penalty, best_pose, best_coords


class SweepConductor:
    """
    High-level orchestrator for rotamer sweeps. Connects pre-docked payloads via
    virtual graph edges to ensure proper rigid-body propagation during collision checks.
    """

    def __init__(self, mol, graph, target_resid, config, chain_id):
        """
        Initializes the SweepConductor.

        :param mol: The combined molecule object.
        :type mol: Mol
        :param graph: The molecular connectivity graph.
        :type graph: MolGraph
        :param target_resid: The target residue sequence number.
        :type target_resid: int
        :param config: The configuration dictionary.
        :type config: dict
        :param chain_id: The target chain identifier.
        :type chain_id: str
        """

        self.mol = mol
        self.graph = graph
        self.resid = target_resid
        self.config = config
        self.chain = chain_id.strip()

        self.prot_sweeper = RotamerSweeper(
            mol, graph, target_resid, config, chain=self.chain
        )

    def optimize(self, moving_atoms, static_atoms):
        """
        Identifies payload molecules via topological bridging, applies the optimal
        steric sweep, and updates the molecular coordinates in place.

        :param moving_atoms: Initial list of moving atoms (typically the sidechain).
        :type moving_atoms: list
        :param static_atoms: Initial list of static environment atoms.
        :type static_atoms: list
        :return: True if a valid pose was calculated and applied, False otherwise.
        :rtype: bool
        """

        scr_atoms = [a for a in self.mol.records if a.res_name.strip() == "LIG"]
        glycan_anchor_idx = None
        scr_anchor_idx = None

        if moving_atoms:
            obj_to_idx = {id(a): i for i, a in enumerate(self.mol.records)}
            hinge_idx = obj_to_idx[id(moving_atoms[0])]
            covalent_tree = nx.node_connected_component(self.graph.nx_graph, hinge_idx)

            patch_atoms = []
            for idx in covalent_tree:
                atom = self.mol.records[idx]
                if (
                    atom.res_name.strip() not in STANDARD_AMINO_ACIDS
                    and atom.res_name.strip() != "LIG"
                ):
                    if atom not in moving_atoms:
                        patch_atoms.append(atom)

            moving_atoms.extend(patch_atoms)

            if scr_atoms and patch_atoms:
                glycan_anchor_idx = obj_to_idx[id(patch_atoms[0])]
                scr_anchor_idx = obj_to_idx[id(scr_atoms[0])]
                self.graph.nx_graph.add_edge(glycan_anchor_idx, scr_anchor_idx)

            if scr_atoms:
                moving_atoms.extend(scr_atoms)

            moving_ids = {id(m) for m in moving_atoms}
            static_atoms = [a for a in static_atoms if id(a) not in moving_ids]

        checker = StericChecker(self.graph, moving_atoms, static_atoms)

        print("\nStarting Sweep...")
        penalty, pose, coords = self.prot_sweeper.get_best_pose(checker, moving_atoms)

        # Remove temporary virtual bridge
        if glycan_anchor_idx is not None and scr_anchor_idx is not None:
            self.graph.nx_graph.remove_edge(glycan_anchor_idx, scr_anchor_idx)

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
