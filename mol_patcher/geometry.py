import numpy as np
from scipy.spatial.transform import Rotation
from typing import List, Tuple
from .mol_record import PdbRecord, Mol
from .utilities import get_distance, get_dihedral
import networkx as nx
import math
import copy
import sys
import json


def get_coords_array(atoms: List[PdbRecord]) -> np.ndarray:
    """
    Extracts an Nx3 numpy array of coordinates from any list of PdbRecord objects.
    """

    return np.array([[a.x, a.y, a.z] for a in atoms])


def rotate_dihedral(atoms, pivot_coord, axis_coord, angle_degrees):
    """
    Applies a quaternion rotation to a list of atoms around a specific bond vector.

    :param atoms: (list) PdbRecord objects to be rotated.
    :param pivot_coord: (tuple/list) The (x,y,z) coordinates of the stationary atom.
    :param axis_coord: (tuple/list) The (x,y,z) coordinates of the atom defining the rotation axis.
    :param angle_degrees: (float) The amount to rotate in degrees.
    :return: (numpy.ndarray) The final transformed coordinate array.
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

    final_pos = None
    for atom in atoms:
        pos = np.array([atom.x, atom.y, atom.z])
        final_pos = rot_obj.apply(pos - pivot) + pivot
        atom.x, atom.y, atom.z = final_pos

    return final_pos


def align_centroids(mol_1, mol_2):
    """
    Calculates the geometric centroid of two atom lists and aligns them.
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
    - If pivot_coord is the glycan centroid, the SCR orbits the glycan.
    - If pivot_coord is the SCR centroid, the SCR spins in place.
    """
    # Normalize the axis of rotation (e.g., [1, 0, 0] for the X-axis)
    axis_arr = np.array(axis_vec)
    axis_norm = axis_arr / np.linalg.norm(axis_arr)

    # Create the rotation object using the standard scipy Rotation matrix
    rot_obj = Rotation.from_rotvec(axis_norm * np.radians(angle_degrees))
    pivot = np.array(pivot_coord)

    # Apply the rotation by entering on the pivot, rotating, and translating back
    for atom in atoms:
        pos = np.array([atom.x, atom.y, atom.z])
        final_pos = rot_obj.apply(pos - pivot) + pivot
        atom.x, atom.y, atom.z = final_pos

    return atoms


def random_rot(atoms, target_centroid, base_coords, max_dist_nm):
    """
    Applies a random 3D rotation + translations to a list of atoms
    """

    # Generate random translation vector within specified max distance
    offset = (np.random.rand(3) - 0.5) * (max_dist_nm * 10) * 2
    target_pos = target_centroid + offset
    # Generate rotation matrix
    rot = Rotation.random().as_matrix()

    # Apply transformation to undistorted base coordinates
    base_centroid = np.mean(base_coords, axis=0)

    temp_coords = []
    for i, atom in enumerate(atoms):
        pos = base_coords[i]
        new_pos = rot @ (pos - base_centroid) + target_pos
        atom.x, atom.y, atom.z = new_pos
        temp_coords.append(new_pos)

    return temp_coords


class PatchAligner:
    """
    Takes a list of equivalent PdbRecord objects from a patch molecule (pfp) and base molecule (protein) as anchors
    Generates the rotational and transformational matrices required to align the anchors
    Apply those operations to all atoms of patch molecule
    Update the coordinates of the patch molecules PdbRecord Objects
    """

    def __init__(
        self,
        patch_atoms: List[PdbRecord],
        patch_anchors: List[PdbRecord],
        target_anchors: List[PdbRecord],
    ):
        """
        Constructs the PatchAligner.

        :param patch_atoms: (list) All PdbRecord objects comprising the patch molecule.
        :param patch_anchors: (list) PdbRecord objects of the patch atoms to be aligned.
        :param target_anchors: (list) PdbRecord objects of the target protein atoms.
        """

        self.patch_atoms = patch_atoms

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

        :param atoms: (list) PdbRecord object(s) comprising the atoms involved in the new bond.
        :param at1: (tuple/list) PdbRecord object for atom on the first molecule that will be replaced.
        :param at2: (tuple/list) PdbRecord object for atom on the first molecule that will be be bonded to at3.
        :param at3: (tuple/list) PdbRecord object for atom on the second molecule that will be be bonded to at2.
        :param at4: (tuple/list) PdbRecord object for atom on the second molecule that will be replaced.
        :param target_bond_length: (float) The length (in angstroms) of the new bond being formed.
        """
        # at1 through at4 are the atoms involved in bonding:
        # initial bonds are at1-at2 and at3-at4 we want to create a bond between at2 and at3 s.t. the coords of at1=at3 and coords of at2=at4

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

        :param p1: (tuple/list) Base primary atom (e.g., Aspargine CA)
        :param p2: (tuple/list) Base anchor atom (e.g., Aspargine N) - Axis of rotation
        :param p3: (tuple/list) Patch anchor atom (e.g., Glycan C1) - Pivot for rotation
        :param p4: (tuple/list) Patch secondary atom (e.g., Glycan Ring O)
        :return: (numpy.ndarray) The final transformed coordinate array.
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

        :param atoms: (list) PdbRecord objects to be rotated (the patch).
        :param p1: (PdbRecord) Base primary atom (e.g., C6)
        :param p2: (PdbRecord) Base anchor atom / The Vertex (e.g., O6)
        :param p3: (PdbRecord) Patch anchor atom (e.g., Galactose C1)
        :param p_ref: (PdbRecord) Reference atom to define the plane (e.g., C5)
        :param target_angle: (float) The desired bond angle in degrees (e.g., 114.0)
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

        # Apply rotation
        for atom in atoms:
            pos = np.array([atom.x, atom.y, atom.z])
            final_pos = best_rot.apply(pos - pivot) + pivot
            atom.x, atom.y, atom.z = final_pos

        return atoms

    def align_from_template(self, base_records, patch_records, template, base_res_idx):
        """
        Executes alignment based on JSON template.
        """
        # Find the connecting residue of the glycan (lowest res_seq)
        root_idx = min(r.res_seq for r in patch_records)

        # Configs
        b_config = template["anchors"]["base"]
        p_config = template["anchors"]["patch"]
        geom = template["geometry"]
        base_res = b_config["res_name"]
        patch_res = p_config["res_name"]

        def get_atom(recs, name, res, seq):
            return next(
                (
                    r
                    for r in recs
                    if r.name.strip() == name
                    and r.res_name.strip() == res
                    and r.res_seq == seq
                ),
                None,
            )

        at_cg = get_atom(base_records, b_config["primary"], base_res, base_res_idx)
        at_nd2 = get_atom(base_records, b_config["anchor"], base_res, base_res_idx)
        at_cb = get_atom(base_records, b_config["ref"], base_res, base_res_idx)
        at_od1 = get_atom(base_records, b_config["carbonyl"], base_res, base_res_idx)

        at_c1 = get_atom(patch_records, p_config["anchor"], patch_res, root_idx)
        at_o5 = get_atom(patch_records, p_config["secondary"], patch_res, root_idx)
        at_o1 = get_atom(patch_records, p_config["leaving"], patch_res, root_idx)

        print(
            f"DEBUG: Patch parsed as -> ResName: '{patch_records[0].res_name}', Seq: {patch_records[0].res_seq}"
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
        self.set_junction_dihedral(
            atoms=patch_records,
            p1=at_od1,
            p2=at_cg,
            p3=at_nd2,
            p4=at_c1,
            target_dih=geom["omega_target"],
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
            p1=at_cg,
            p2=at_nd2,
            p3=at_c1,
            p4=at_o5,
            target_dih=geom["phi_target"],
        )

        return patch_records


class MolGraph:
    def __init__(self, mol: Mol, distXX=1.65, distXH=1.15, itp_to_pdb_map=None):
        """
        Creates a connectivity matrix of the molecule. A connectivity matrix holds the information of which atoms are bonded and to what.

        :param distXX: The max distance between two atoms (not hydrogen) to be considered a bond
        :param distXH: The max distance between any atom and a hydrogen atom to be considered a bond
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
        """
        new_instance = copy.copy(self)
        if self.nx_graph is not None:
            new_instance.nx_graph = self.nx_graph.copy()
        return new_instance

    def _build_edges_from_topology(self, itp_to_pdb_map):
        """
        Builds the connectivity matrix and graph directly from the ITP bonds list,
        bypassing spatial distance calculations entirely.
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
        Executes distance checks and builds the connectivity matrix.

        To avoid protein-wide O(N^2) scaling and inter-residue clashes, distance
        calculations are localized to atoms within the same residue. The peptide
        backbone is then stitched together using an O(1) dictionary lookup map.

        :return: None. Updates the internal self.matrix and self.nx_graph variables.
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
    def __init__(self, mol_graph, moving_atoms, static_atoms, tolerance=0.55):
        self.mol_graph = mol_graph

        # Filter out atoms that do not move during sidechain sweeps.
        # This prevents permanent unresolvable native clashes from paralyzing the script.
        self.moving_atoms = moving_atoms
        self.static_atoms = static_atoms
        self.tolerance = tolerance

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
                self.mol_graph.nx_graph, m_idx, cutoff=12
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

    def check_clashes(self, limit_to_atoms=None):
        """
        Evaluates distances between moving and static atoms.
        :param limit_to_atoms: (list, optional) If provided, only check clashes originating from these PdbRecord objects.
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

                m2_element = m_atom_2.name.strip()
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
        Calculates a penalty score for the current pose by summing the penetration depth of all overlaps.
        """
        check_list = limit_to_atoms if limit_to_atoms is not None else self.moving_atoms
        environment_atoms = self.static_atoms.copy()
        if limit_to_atoms is not None:
            environment_atoms.extend(
                [a for a in self.moving_atoms if a not in limit_to_atoms]
            )

        total_penalty = 0.0
        clash_count = 0

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
                    total_penalty += threshold - actual_dist
                    clash_count += 1

        return total_penalty, clash_count
