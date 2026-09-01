"""
The spatial and topological engine for MolPatcher.
Handles 3D coordinate manipulation, matrix alignments, graph building, and steric evaluations.
"""

import numpy as np
from scipy.spatial.transform import Rotation
from typing import List
from .mol_record import PdbRecord, Mol
from .utilities import get_distance, get_dihedral
import networkx as nx
import copy
import sys


def get_coords_array(atoms: List[PdbRecord]) -> np.ndarray:
    """
    Extracts an Nx3 numpy array of coordinates from any list of PdbRecord objects.

    :param atoms: List of atoms to extract coordinates from.
    :type atoms: list
    :return: An Nx3 matrix of (x, y, z) coordinates.
    :rtype: numpy.ndarray
    """

    return np.array([[a.x, a.y, a.z] for a in atoms])


def rotate_dihedral(atoms, pivot_coord, axis_coord, angle_degrees):
    """
    Applies a quaternion rotation to a list of atoms around a specific bond vector.

    :param atoms: Objects to be rotated.
    :type atoms: list
    :param pivot_coord: The (x, y, z) coordinates of the stationary atom.
    :type pivot_coord: tuple or list
    :param axis_coord: The (x, y, z) coordinates of the atom defining the rotation axis.
    :type axis_coord: tuple or list
    :param angle_degrees: The amount to rotate in degrees.
    :type angle_degrees: float
    :raises ValueError: If pivot and axis coordinates are identical.
    :return: The final transformed coordinate array.
    :rtype: numpy.ndarray
    """

    # Define and normalize the axis vector
    axis_vec = np.array(axis_coord) - np.array(pivot_coord)
    axis_norm = np.linalg.norm(axis_vec)

    if axis_norm == 0:
        raise ValueError(
            "Pivot and axis coordinates are identical; cannot define rotation axis."
        )

    axis_vec /= axis_norm

    # Create the rotation object using the rotation vector (angle * axis)
    rot_obj = Rotation.from_rotvec(axis_vec * np.radians(angle_degrees))
    pivot = np.array(pivot_coord)

    # Transform each atom
    if not atoms:
        return None  # Return None if there is nothing to rotate

    # Vectorized Rotation
    coords = get_coords_array(atoms)
    final_coords = rot_obj.apply(coords - pivot) + pivot

    return final_coords


def align_centroids(mol_1, mol_2):
    """
    Calculates the geometric centroid of two atom lists and aligns them.

    :param mol_1: The reference molecule atoms.
    :type mol_1: list
    :param mol_2: The moving molecule atoms to translate.
    :type mol_2: list
    :return: The translated list of moving atoms.
    :rtype: list
    """

    mol_1_coords = np.array([[a.x, a.y, a.z] for a in mol_1])
    mol_2_coords = np.array([[a.x, a.y, a.z] for a in mol_2])

    mol_1_centroid = np.mean(mol_1_coords, axis=0)
    mol_2_centroid = np.mean(mol_2_coords, axis=0)

    translation_vector = mol_1_centroid - mol_2_centroid

    for atom in mol_2:
        atom.x += translation_vector[0]
        atom.y += translation_vector[1]
        atom.z += translation_vector[2]

    return mol_2


def translate_along_vector(atoms, start_coord, end_coord, distance=1.5):
    """
    Pushes a list of atoms a specific distance along a defined directional vector.

    :param atoms: Objects to be translated.
    :type atoms: list
    :param start_coord: The starting (x, y, z) coordinates.
    :type start_coord: tuple or list
    :param end_coord: The ending (x, y, z) coordinates.
    :type end_coord: tuple or list
    :param distance: The translation distance in Angstroms, defaults to 1.5.
    :type distance: float, optional
    :raises ValueError: If start and end coordinates are identical.
    :return: The translated list of atoms.
    :rtype: list
    """

    direction_vec = np.array(end_coord) - np.array(start_coord)
    vec_norm = np.linalg.norm(direction_vec)

    if vec_norm == 0:
        raise ValueError("Start and end coordinates are identical.")

    # Normalize vector to length 1, multiply by target distance (1.5 Å)
    normalized_vec = direction_vec / vec_norm
    translation = normalized_vec * distance

    for atom in atoms:
        atom.x += translation[0]
        atom.y += translation[1]
        atom.z += translation[2]

    return atoms


def rotate_around_pivot(atoms, pivot_coord, axis_vec, angle_degrees):
    """
    Rotates a list of atoms around a specific pivot point in 3D space.

    :param atoms: Objects to be rotated.
    :type atoms: list
    :param pivot_coord: The (x, y, z) coordinates of the rotation center.
    :type pivot_coord: tuple or list
    :param axis_vec: The (x, y, z) vector defining the axis of rotation.
    :type axis_vec: tuple or list
    :param angle_degrees: The rotation amount in degrees.
    :type angle_degrees: float
    :return: The rotated list of atoms.
    :rtype: list
    """

    # Normalize the axis of rotation (e.g., [1, 0, 0] for the X-axis)
    axis_arr = np.array(axis_vec)
    axis_norm = axis_arr / np.linalg.norm(axis_arr)

    # Create the rotation object using the standard scipy Rotation matrix
    rot_obj = Rotation.from_rotvec(axis_norm * np.radians(angle_degrees))
    pivot = np.array(pivot_coord)

    # Vectorized rotation
    coords = get_coords_array(atoms)
    final_coords = rot_obj.apply(coords - pivot) + pivot

    for i, atom in enumerate(atoms):
        atom.x, atom.y, atom.z = final_coords[i]

    return atoms


class PatchAligner:
    """
    Generates rotational and transformational matrices to align equivalent anchors
    between a patch molecule and a base protein.
    """

    def __init__(
        self,
        patch_atoms: List[PdbRecord],
        patch_anchors: List[PdbRecord],
        target_anchors: List[PdbRecord],
    ):
        """
        Constructs the PatchAligner.

        :param patch_atoms: All PdbRecord objects comprising the patch molecule.
        :type patch_atoms: list
        :param patch_anchors: PdbRecord objects of the patch atoms to be aligned.
        :type patch_anchors: list
        :param target_anchors: PdbRecord objects of the target protein atoms.
        :type target_anchors: list
        """

        self.patch_atoms = patch_atoms
        self.target_chain = target_anchors[0].chain if target_anchors else " "

        if not patch_anchors or not target_anchors:
            self.rotation_object = None
            self.translation_vector = None
            return

        # Extract coords from provided atom objects
        patch_anchor_coords = get_coords_array(patch_anchors)
        target_anchor_coords = get_coords_array(target_anchors)

        patch_centroid = np.mean(patch_anchor_coords, axis=0)
        target_centroid = np.mean(target_anchor_coords, axis=0)

        patch_centered = patch_anchor_coords - patch_centroid
        target_centered = target_anchor_coords - target_centroid

        self.rotation_object, rmsd, *_ = Rotation.align_vectors(
            target_centered, patch_centered
        )
        rotated_patch_centroid = self.rotation_object.apply(patch_centroid)
        self.translation_vector = target_centroid - rotated_patch_centroid

    def implement_align(self):
        """
        Applies the calculated alignment matrices to the entire patch molecule.

        :return: The transformed list of patch atoms.
        :rtype: list
        """

        if self.rotation_object is None or self.translation_vector is None:
            return self.patch_atoms

        coords = np.array([[a.x, a.y, a.z] for a in self.patch_atoms])
        rotated_coords = self.rotation_object.apply(coords)
        final_coords = rotated_coords + self.translation_vector

        for i, atom in enumerate(self.patch_atoms):
            atom.x, atom.y, atom.z = final_coords[i]

        return self.patch_atoms

    def align_single_bond(
        self, atoms: List[PdbRecord], at1, at2, at3, at4, target_bond_length
    ):
        """
        Translates and rotates the patch molecule using a single bond vector alignment.

        :param atoms: PdbRecord objects comprising the patch molecule.
        :type atoms: list
        :param at1: Base primary atom that will be replaced.
        :type at1: PdbRecord
        :param at2: Base anchor atom forming the bond.
        :type at2: PdbRecord
        :param at3: Patch anchor atom forming the bond.
        :type at3: PdbRecord
        :param at4: Patch leaving group atom that will be replaced.
        :type at4: PdbRecord
        :param target_bond_length: The length of the new bond in Angstroms.
        :type target_bond_length: float
        :return: The transformed list of patch atoms.
        :rtype: list
        """

        # at1 through at4 are the atoms involved in bonding:
        # initial bonds are at1-at2 and at3-at4 we want to create a bond between at2 and at3 such that the coords of at1=at3 and coords of at2=at4

        self.at1_coords = np.array([at1.x, at1.y, at1.z])
        self.at2_coords = np.array([at2.x, at2.y, at2.z])
        self.at3_coords = np.array([at3.x, at3.y, at3.z])
        self.at4_coords = np.array([at4.x, at4.y, at4.z])

        bond1_vec = self.at2_coords - self.at1_coords
        bond2_vec = self.at4_coords - self.at3_coords

        target_vec = bond1_vec
        target_magnitude = np.linalg.norm(target_vec)
        target_vec_normalized = target_vec / target_magnitude

        target_at3_pos = self.at2_coords + (target_vec_normalized * target_bond_length)

        bond_rotation_obj, rmsd, *_ = Rotation.align_vectors([-target_vec], [bond2_vec])

        for atom in atoms:
            coords = np.array([atom.x, atom.y, atom.z])
            centered_coords = coords - self.at3_coords
            rotated_coords = bond_rotation_obj.apply(centered_coords)
            translated_coords = rotated_coords + target_at3_pos
            atom.x, atom.y, atom.z = translated_coords
        return atoms

    def set_junction_dihedral(self, atoms: List[PdbRecord], p1, p2, p3, p4, target_dih):
        """
        Spins the patch molecule around the new junction bond to match a specific dihedral angle.

        :param atoms: The patch atoms to be rotated.
        :type atoms: list
        :param p1: Base primary atom.
        :type p1: PdbRecord
        :param p2: Base anchor atom (axis of rotation).
        :type p2: PdbRecord
        :param p3: Patch anchor atom (pivot for rotation).
        :type p3: PdbRecord
        :param p4: Patch secondary atom.
        :type p4: PdbRecord
        :param target_dih: The target dihedral angle in degrees.
        :type target_dih: float
        :return: The transformed list of patch atoms.
        :rtype: list
        """

        self.p1_coords = np.array([p1.x, p1.y, p1.z])
        self.p2_coords = np.array([p2.x, p2.y, p2.z])
        self.p3_coords = np.array([p3.x, p3.y, p3.z])
        self.p4_coords = np.array([p4.x, p4.y, p4.z])

        current_dih = get_dihedral(
            self.p1_coords, self.p2_coords, self.p3_coords, self.p4_coords
        )
        delta = target_dih - current_dih
        rotate_dihedral(
            atoms, self.p2_coords, self.p3_coords, delta
        )  # fixed so it's properly p2 -> p3 axis

        return atoms

    def set_junction_angle(
        self,
        atoms: List[PdbRecord],
        p1,
        p2,
        p3,
        p_ref,
        target_angle,
        repulsion_coord=None,
    ):
        """
        Bends the bond angle of a junction using a reference atom to define the bending plane.

        :param atoms: The patch atoms to be bent.
        :type atoms: list
        :param p1: Base primary atom.
        :type p1: PdbRecord
        :param p2: Base anchor atom / The vertex.
        :type p2: PdbRecord
        :param p3: Patch anchor atom.
        :type p3: PdbRecord
        :param p_ref: Reference atom to define the bending plane.
        :type p_ref: PdbRecord
        :param target_angle: The desired bond angle in degrees.
        :type target_angle: float
        :param repulsion_coord: Optional coordinate to maximize distance from during bending.
        :type repulsion_coord: list or tuple, optional
        :raises ValueError: If reference atoms are collinear.
        :return: The transformed list of patch atoms.
        :rtype: list
        """

        v1 = np.array([p1.x, p1.y, p1.z]) - np.array([p2.x, p2.y, p2.z])
        v3 = np.array([p3.x, p3.y, p3.z]) - np.array([p2.x, p2.y, p2.z])
        v_ref = np.array([p_ref.x, p_ref.y, p_ref.z]) - np.array([p2.x, p2.y, p2.z])

        cosine_angle = np.dot(v1, v3) / (np.linalg.norm(v1) * np.linalg.norm(v3))
        cosine_angle = np.clip(cosine_angle, -1.0, 1.0)
        current_angle = np.degrees(np.arccos(cosine_angle))

        delta = target_angle - current_angle
        hinge_axis = np.cross(v1, v_ref)
        axis_norm = np.linalg.norm(hinge_axis)

        if axis_norm == 0:
            raise ValueError(
                "Reference atoms are collinear; cannot define a bending plane."
            )

        hinge_axis_normalized = hinge_axis / axis_norm
        pivot = np.array([p2.x, p2.y, p2.z])

        rot_forward = Rotation.from_rotvec(hinge_axis_normalized * np.radians(delta))

        if repulsion_coord is not None:
            # Calc the opposite rotation
            rot_backward = Rotation.from_rotvec(
                hinge_axis_normalized * np.radians(-delta)
            )

            # Find the patch center of mass
            patch_coords = np.array([[a.x, a.y, a.z] for a in atoms])
            patch_centroid = np.mean(patch_coords, axis=0)

            # Test centroid under both hinge rotations
            test_pos_1 = rot_forward.apply(patch_centroid - pivot) + pivot
            test_pos_2 = rot_backward.apply(patch_centroid - pivot) + pivot

            # Measure distances between the new patch centroids and the repulsion coordinate
            dist1 = np.linalg.norm(test_pos_1 - np.array(repulsion_coord))
            dist2 = np.linalg.norm(test_pos_2 - np.array(repulsion_coord))

            # Select the rotation that maximizes the distance
            best_rot = rot_forward if dist1 > dist2 else rot_backward
        else:
            best_rot = rot_forward

        # Vectorized bending
        coords = get_coords_array(atoms)
        final_coords = best_rot.apply(coords - pivot) + pivot

        for i, atom in enumerate(atoms):
            atom.x, atom.y, atom.z = final_coords[i]

        return atoms

    def align_from_template(self, base_records, patch_records, template, base_res_idx):
        """
        Executes sequence alignment based on predefined values in the JSON template.

        :param base_records: PdbRecord objects of the base protein.
        :type base_records: list
        :param patch_records: PdbRecord objects of the patch molecule.
        :type patch_records: list
        :param template: The configuration dictionary.
        :type template: dict
        :param base_res_idx: The sequence index of the target residue.
        :type base_res_idx: int
        :return: The transformed list of patch atoms.
        :rtype: list
        """

        root_idx = min(r.res_seq for r in patch_records)

        # Configs
        b_config = template["anchors"]["base"]
        p_config = template["anchors"]["patch"]
        geom = template["geometry"]
        base_res = b_config["res_name"]
        patch_res = p_config["res_name"]

        def get_atom(recs, name, res, seq, chain):
            return next(
                (
                    r
                    for r in recs
                    if r.name.strip() == name
                    and r.res_name.strip() == res
                    and r.res_seq == seq
                    and r.chain.strip() == chain.strip()
                ),
                None,
            )

        at_cg = get_atom(
            base_records, b_config["primary"], base_res, base_res_idx, self.target_chain
        )
        at_nd2 = get_atom(
            base_records, b_config["anchor"], base_res, base_res_idx, self.target_chain
        )
        at_cb = get_atom(
            base_records, b_config["ref"], base_res, base_res_idx, self.target_chain
        )

        patch_chain = patch_records[0].chain

        at_c1 = get_atom(
            patch_records, p_config["anchor"], patch_res, root_idx, patch_chain
        )
        at_o5 = get_atom(
            patch_records, p_config["secondary"], patch_res, root_idx, patch_chain
        )
        at_o1 = get_atom(
            patch_records, p_config["leaving"], patch_res, root_idx, patch_chain
        )

        print(
            f"Patch parsed as -> ResName: '{patch_records[0].res_name}', Seq: {patch_records[0].res_seq}"
        )

        # Execute the geometric alignment
        self.align_single_bond(
            atoms=patch_records,
            at1=at_cg,
            at2=at_nd2,
            at3=at_c1,
            at4=at_o1,
            target_bond_length=geom["bond_length"],
        )
        self.set_junction_angle(
            atoms=patch_records,
            p1=at_cg,
            p2=at_nd2,
            p3=at_c1,
            p_ref=at_cb,
            target_angle=geom["angle_target"],
        )

        self.set_junction_dihedral(
            atoms=patch_records,
            p1=at_cb,
            p2=at_cg,
            p3=at_nd2,
            p4=at_c1,
            target_dih=geom["omega_target"],
        )
        self.set_junction_dihedral(
            atoms=patch_records,
            p1=at_cg,
            p2=at_nd2,
            p3=at_c1,
            p4=at_o5,
            target_dih=geom["phi_target"],
        )

        return patch_records


class MolGraph:
    """
    Constructs a NetworkX connectivity matrix of the molecule for topology traversals.
    """

    def __init__(self, mol: Mol, distXX=1.65, distXH=1.15, itp_to_pdb_map=None):
        """
        Initializes the connectivity matrix mapping bonded interactions.

        :param mol: The standardized molecule container.
        :type mol: Mol
        :param distXX: The maximum distance between two heavy atoms to be considered bonded.
        :type distXX: float
        :param distXH: The maximum distance between a heavy atom and hydrogen to be considered bonded.
        :type distXH: float
        :param itp_to_pdb_map: Optional mapping dictionary to bypass spatial calculation.
        :type itp_to_pdb_map: dict, optional
        """

        self.mol = mol
        self.distXX = distXX
        self.distXH = distXH

        # State variables
        self.matrix = None
        self.nx_graph = nx.Graph()
        self.Nmols = None
        self.is_distance_based = False

        for i, record in enumerate(mol.records):
            self.nx_graph.add_node(i)
        if itp_to_pdb_map and mol.bonds:
            self._build_edges_from_topology(itp_to_pdb_map)
        else:
            self._build_matrix()  # Automatically build the graph upon initialization
            self.is_distance_based = True

    def copy(self):
        """
        Creates a shallow copy of the MolGraph.

        :return: A copied instance of the graph.
        :rtype: MolGraph
        """
        new_instance = copy.copy(self)
        if self.nx_graph is not None:
            new_instance.nx_graph = self.nx_graph.copy()
        return new_instance

    def _build_edges_from_topology(self, itp_to_pdb_map):
        """
        Builds the graph directly from the ITP bonds list, bypassing spatial distance sweeps.
        """

        Natoms = len(self.mol.records)
        conn_mat = np.zeros((Natoms, Natoms))

        # Create dictionary mapping the PdbRecord's memory ID to its list index
        record_to_idx = {id(r): i for i, r in enumerate(self.mol.records)}

        for bond in self.mol.bonds:
            # Look up PdbRecord objects using ITP topology numbers
            rec1 = itp_to_pdb_map.get(bond.a1)
            rec2 = itp_to_pdb_map.get(bond.a2)

            # If both atoms exist in the PDB (weren't deleted during patching)
            if rec1 and rec2:
                # Find exact indices in the mol.records list
                idx1 = record_to_idx.get(id(rec1))
                idx2 = record_to_idx.get(id(rec2))

                if idx1 is not None and idx2 is not None:
                    # Add the connection to both the NetworkX graph and numpy matrix
                    self.nx_graph.add_edge(idx1, idx2)
                    conn_mat[idx1, idx2] = 1
                    conn_mat[idx2, idx1] = 1

        # Store final state variables
        self.matrix = conn_mat
        self.Nmols = nx.number_connected_components(self.nx_graph)

    def _build_matrix(self):
        """
        Executes localized distance checks to map the backbone and build the connectivity matrix.
        """

        Natoms = len(self.mol.records)
        conn_mat = np.zeros((Natoms, Natoms))

        def get_coords(idx):
            return [
                self.mol.records[idx].x,
                self.mol.records[idx].y,
                self.mol.records[idx].z,
            ]

        def is_hydrogen(idx):
            return self.mol.records[idx].name.strip().startswith("H")

        # Creates O(1) lookup dictionary for the exact backbone Carbonyl (C) and Amide (N) atoms.
        # Allows the peptide chain to be manually linked later without global distance sweeps.
        backbone_map = {}

        spinner = ["|", "/", "-", "\\"]
        for i in range(Natoms):
            # --- SPINNER (Indicates code is running)---
            if i % 30 == 0 or i == Natoms - 1:
                char = spinner[(i // 30) % 4]
                sys.stdout.write(f"\rBuilding connectivity matrix... {char}")
                sys.stdout.flush()
            # ---------------------------------

            # --- MAP THE PROTEIN BACKBONE ---
            rec_i = self.mol.records[i]

            atom_name = rec_i.name.strip()
            if atom_name in ["C", "N"]:  # Check if atom is a connecting carbonyl/amide
                key = (
                    rec_i.chain,
                    rec_i.res_seq,
                )  # If it is, save it to the dictionary
                if key not in backbone_map:
                    backbone_map[key] = {}
                backbone_map[key][atom_name] = i

            for j in range(
                i + 1, Natoms
            ):  # i + 1 prevents self-checking and double-counting
                rec_j = self.mol.records[j]

                if rec_i.res_seq == rec_j.res_seq and rec_i.chain == rec_j.chain:
                    dist = get_distance(get_coords(i), get_coords(j))

                    is_h = is_hydrogen(i) or is_hydrogen(j)

                    if is_h and dist <= self.distXH:
                        conn_mat[i, j] = 1
                        conn_mat[j, i] = 1

                    elif not is_h and dist <= self.distXX:
                        conn_mat[i, j] = 1
                        conn_mat[j, i] = 1

        sys.stdout.write("\b \n")
        sys.stdout.flush()

        # --- CONNECT THE PROTEIN BACKBONE ---
        # Iterate through mapped residues to connect C of res n to N of res n+1
        for (chain, res_seq), atoms in backbone_map.items():
            next_key = (chain, res_seq + 1)

            if next_key in backbone_map:
                c_idx = atoms.get("C")
                n_idx = backbone_map[next_key].get("N")

                # If both atoms exist, generate a bond in the matrix
                if c_idx is not None and n_idx is not None:
                    # Double check that the atoms are close enough to actually bond
                    if get_distance(get_coords(c_idx), get_coords(n_idx)) < 2.0:
                        conn_mat[c_idx, n_idx] = 1
                        conn_mat[n_idx, c_idx] = 1

        # Remove bifurcated H atoms
        for i in range(Natoms):
            if is_hydrogen(i) and np.sum(conn_mat[i, :]) > 1:
                bonded_indices = np.where(conn_mat[i, :] == 1)[0]

                heavy_bonded = [idx for idx in bonded_indices if not is_hydrogen(idx)]

                # Fallback if it is only bonded to hydrogens
                if not heavy_bonded:
                    heavy_bonded = bonded_indices

                closest_idx = None
                min_dist = float("inf")

                # Find the closest atom
                for b_idx in heavy_bonded:
                    dist = get_distance(get_coords(i), get_coords(b_idx))
                    if dist < min_dist:
                        min_dist = dist
                        closest_idx = b_idx

                # Wipe the hydrogen's connections, then restore only the closest one
                if closest_idx is not None:
                    conn_mat[i, :] = 0
                    conn_mat[:, i] = 0
                    conn_mat[i, closest_idx] = 1
                    conn_mat[closest_idx, i] = 1

        # Store final state variables
        self.matrix = conn_mat
        self.nx_graph = nx.Graph(conn_mat)
        self.Nmols = nx.number_connected_components(self.nx_graph)


class StericChecker:
    """
    Evaluates spatial environments using a localized graph-traversal search to identify overlaps based on VdW radii.
    """

    def __init__(
        self, mol_graph, moving_atoms, static_atoms, tolerance=0.55, docking_target=None
    ):
        """
        Initializes the StericChecker.

        :param mol_graph: The molecular connectivity graph.
        :type mol_graph: MolGraph
        :param moving_atoms: Atoms that will be transformed during sweeps.
        :type moving_atoms: list
        :param static_atoms: The rigid environment atoms.
        :type static_atoms: list
        :param tolerance: Fraction of combined VdW radii allowed before triggering a clash.
        :type tolerance: float
        :param docking_target: Optional target coordinates for localized spatial filtering.
        :type docking_target: list, optional
        """

        self.mol_graph = mol_graph

        # Filter out atoms that do not move during sidechain sweeps.
        # This prevents permanent unresolvable native clashes from paralyzing the script.
        self.moving_atoms = moving_atoms
        self.static_atoms = static_atoms
        self.tolerance = tolerance

        # Filter for non-covalent interaction tolerance
        if docking_target is not None and len(docking_target) > 0:
            target_coords = np.array([[a.x, a.y, a.z] for a in docking_target])
            centroid = np.mean(target_coords, axis=0)

            self.static_atoms = [
                s
                for s in static_atoms
                if get_distance((s.x, s.y, s.z), centroid) < 30.0
            ]
            print(
                f"    [Spatial Filter] Reduced search space from {len(static_atoms)} to {len(self.static_atoms)} atoms."
            )

        else:
            self.static_atoms = static_atoms

        # Effective/bonded van der Waals radii (in Angstroms). See J. Charry, A. Tkatchenko, J. Chem. Theory Comput. 2024, 20, 7469–7478.
        self.vdw_radii = {
            "H": 1.10,  # modified from 1.20 since H's are artificially added in the first place
            "C": 1.70,
            "N": 1.55,
            "O": 1.52,
            "P": 1.80,
            "S": 1.80,
            "F": 1.47,  # Bondi Radius. See A. Bondi, J. Phys. Chem. 1964, 68, 441–451.
        }

        self.atom_to_idx = {
            atom.serial: i for i, atom in enumerate(self.mol_graph.mol.records)
        }
        self.excluded_pairs = set()

        for m_atom in self.moving_atoms:
            m_idx = self.atom_to_idx[m_atom.serial]

            close_neighbors = nx.single_source_shortest_path_length(
                self.mol_graph.nx_graph, m_idx, cutoff=3
            )

            for s_atom in self.static_atoms:
                s_idx = self.atom_to_idx[s_atom.serial]
                # O(1) dictionary lookup
                if s_idx in close_neighbors:
                    self.excluded_pairs.add((m_atom.serial, s_atom.serial))

            for m_atom_2 in self.moving_atoms:
                m2_idx = self.atom_to_idx[m_atom_2.serial]
                if m2_idx in close_neighbors:
                    self.excluded_pairs.add((m_atom.serial, m_atom_2.serial))

        self.res_depths = {}
        if self.moving_atoms:
            hinge_idx = self.atom_to_idx[self.moving_atoms[0].serial]
            self.res_depths[hinge_idx] = 0

            for u, v in nx.bfs_edges(self.mol_graph.nx_graph, hinge_idx):
                u_rec = self.mol_graph.mol.records[u]
                v_rec = self.mol_graph.mol.records[v]

                if u_rec.res_seq == v_rec.res_seq:
                    self.res_depths[v] = self.res_depths.get(u, 0)
                else:
                    self.res_depths[v] = self.res_depths.get(u, 0) + 1

    def check_clashes(self, limit_to_atoms=None):
        """
        Evaluates distances between moving and static atoms and returns immediately upon finding a clash.

        :param limit_to_atoms: Specific PdbRecord objects to limit the search scope, defaults to None.
        :type limit_to_atoms: list, optional
        :return: A tuple of the actual clash distance and the formatting error message, or (None, None) if clean.
        :rtype: tuple
        """

        # If no list provided, check all known moving atoms defined during init
        check_list = limit_to_atoms if limit_to_atoms is not None else self.moving_atoms

        # Create a complete environment list that includes static atoms + any sidechain atoms that are currently stationary
        environment_atoms = self.static_atoms.copy()
        if limit_to_atoms is not None:
            environment_atoms.extend(
                [a for a in self.moving_atoms if a not in limit_to_atoms]
            )

        for m_atom in check_list:
            m_element = m_atom.name.strip()[0]
            m_vdw = self.vdw_radii.get(m_element, 1.50)

            for s_atom in environment_atoms:
                if m_atom.serial == s_atom.serial:
                    continue
                if (m_atom.serial, s_atom.serial) in self.excluded_pairs:
                    continue
                if s_atom.name.strip().startswith("H"):
                    continue

                s_element = s_atom.name.strip()[0]
                s_vdw = self.vdw_radii.get(s_element, 1.50)

                actual_dist = get_distance(
                    (m_atom.x, m_atom.y, m_atom.z), (s_atom.x, s_atom.y, s_atom.z)
                )

                threshold = (m_vdw + s_vdw) * self.tolerance
                if actual_dist < threshold:
                    msg = f"CLASH: {m_atom.name:<4} hit {s_atom.name:<4} ({s_atom.res_seq:<3}) | Dist: {actual_dist:.2f} Å"
                    return actual_dist, msg

            for m_atom_2 in check_list:
                if m_atom.serial == m_atom_2.serial:
                    continue
                if (m_atom.serial, m_atom_2.serial) in self.excluded_pairs:
                    continue
                if m_atom_2.name.strip().startswith("H"):
                    continue

                m2_element = m_atom_2.name.strip()[0]
                m2_vdw = self.vdw_radii.get(m2_element, 1.50)

                actual_dist = get_distance(
                    (m_atom.x, m_atom.y, m_atom.z), (m_atom_2.x, m_atom_2.y, m_atom_2.z)
                )

                threshold = (m_vdw + m2_vdw) * self.tolerance
                if actual_dist < threshold:
                    msg = f"CLASH: {m_atom.name:<4} hit {m_atom_2.name:<4} ({m_atom_2.res_seq:<3}) | Dist: {actual_dist:.2f} Å"
                    return actual_dist, msg

        return None, None  # No clash

    def score_pose(self, limit_to_atoms=None):
        """
        Calculates a penalty score using vectorized numpy broadcasting to evaluate penetration depths.

        :param limit_to_atoms: Specific PdbRecord objects to score, defaults to None.
        :type limit_to_atoms: list, optional
        :return: A tuple containing the total accumulated penalty and the number of distinct clashes.
        :rtype: tuple of (float, int)
        """

        check_list = limit_to_atoms if limit_to_atoms is not None else self.moving_atoms

        base_environment = self.static_atoms.copy()
        if limit_to_atoms is not None:
            base_environment.extend(
                [a for a in self.moving_atoms if a not in limit_to_atoms]
            )

        env_heavy = [s for s in base_environment if not s.name.strip().startswith("H")]

        moving_coords = np.array([[a.x, a.y, a.z] for a in check_list])
        env_coords = np.array([[s.x, s.y, s.z] for s in env_heavy])

        m_vdw = np.array(
            [self.vdw_radii.get(a.name.strip()[0], 1.50) for a in check_list]
        )
        e_vdw = np.array(
            [self.vdw_radii.get(s.name.strip()[0], 1.50) for s in env_heavy]
        )

        diff = moving_coords[:, np.newaxis, :] - env_coords[np.newaxis, :, :]
        distances = np.linalg.norm(diff, axis=2)

        thresholds = (m_vdw[:, np.newaxis] + e_vdw[np.newaxis, :]) * self.tolerance

        # Glycodiside penalties decrease we move to each succesive sugar. This accounts for glycosidic bond flexibility.
        depths = np.array(
            [self.res_depths.get(self.atom_to_idx[a.serial], 0) for a in check_list]
        )
        weights = 1.0 / (depths + 1.0)

        penetrations = np.where(distances < thresholds, thresholds - distances, 0.0)

        weighted_penetrations = penetrations * weights[:, np.newaxis]
        penalties = (weighted_penetrations**3) * 50

        for i, m_atom in enumerate(check_list):
            for j, s_atom in enumerate(env_heavy):
                if (
                    m_atom.serial == s_atom.serial
                    or (m_atom.serial, s_atom.serial) in self.excluded_pairs
                ):
                    penalties[i, j] = 0.0

        total_penalty = np.sum(penalties)
        clash_count = np.count_nonzero(penalties)

        # Clash Debugging
        # clash_indices = np.argwhere(penalties > 0)
        # if len(clash_indices) > 0:
        #     print("\n--- DEBUG: Breakdown of Pose Penalty ---")
        #     for i, j in clash_indices:
        #         m_atom = check_list[i]
        #         s_atom = env_heavy[j]
        #         print(
        #             f"  {m_atom.name:<4} ({m_atom.res_seq}) hit {s_atom.name:<4} ({s_atom.res_seq}) | Added Penalty: {penalties[i, j]:.2f}"
        #         )
        #     print("----------------------------------------")

        return total_penalty, clash_count

    def score_docking(self, limit_to_atoms=None, contact_threshold=5.5):
        """
        Calculates a distance score for docking using Numpy arrays.

        :param limit_to_atoms: Specific PdbRecord objects to score, defaults to None.
        :type limit_to_atoms: list, optional
        :param contact_threshold: Max distance in Angstroms to be considered a positive contact.
        :type contact_threshold: float, optional
        :return: A float score where lower (more negative) is better.
        :rtype: float
        """

        check_list = limit_to_atoms if limit_to_atoms is not None else self.moving_atoms

        m_coords = np.array([[a.x, a.y, a.z] for a in check_list])
        s_coords = np.array([[a.x, a.y, a.z] for a in self.static_atoms])

        m_vdw = np.array(
            [self.vdw_radii.get(a.name.strip()[0], 1.5) for a in check_list]
        )
        s_vdw = np.array(
            [self.vdw_radii.get(a.name.strip()[0], 1.5) for a in self.static_atoms]
        )

        # use np.newaxis for broadcasting -> adding dimension(s) of size 1 so that we can perform element operations
        diff = m_coords[:, np.newaxis, :] - s_coords[np.newaxis, :, :]
        distances = np.linalg.norm(
            diff, axis=2
        )  # axis=2 so that we use standard 3D distance formula along 3rd dimension.

        docking_tolerance = 0.9
        clash_thresholds = (
            m_vdw[:, np.newaxis] + s_vdw[np.newaxis, :]
        ) * docking_tolerance

        clash_mask = distances < clash_thresholds
        total_penalty = np.sum(
            (clash_thresholds[clash_mask] - distances[clash_mask]) * 100
        )

        contact_mask = (distances >= clash_thresholds) & (distances < contact_threshold)
        contact_reward = np.sum(contact_mask) * 0.1

        return total_penalty - contact_reward
