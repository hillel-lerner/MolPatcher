import numpy as np
from scipy.spatial.transform import Rotation
from typing import List, Tuple
from .mol_record import PdbRecord, Mol
from .utilities import get_distance, get_dihedral
from .mol_record import Mol
import networkx as nx
import math
import copy
import sys

    
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
            raise ValueError("Pivot and axis coordinates are identical; cannot define rotation axis.")
        
        axis_vec /= axis_norm
        
        # Create the rotation object using the rotation vector (angle * axis)
        rot_obj = Rotation.from_rotvec(axis_vec * np.radians(angle_degrees))  # Scipy's from_rotvec uses quaternions for the calculation
        pivot = np.array(pivot_coord)
        
        # Transform each atom
        if not atoms:
            return None # Return None if there is nothing to rotate
            
        final_pos = None
        for atom in atoms:
            pos = np.array([atom.x, atom.y, atom.z])
            final_pos = rot_obj.apply(pos - pivot) + pivot
            atom.x, atom.y, atom.z = final_pos
            
        return final_pos

class PatchAligner:

    """
    Takes a list of equivalent PdbRecord objects from a patch molecule (pfp) and base molecule (protein) as anchors
    Generates the rotational and transformational matrices required to align the anchors
    Apply those operations to all atoms of patch molecule
    Update the coordinates of the patch molecules PdbRecord Objects
    """

    def __init__(self, pfp_atoms: List[PdbRecord], pfp_anchors: List[PdbRecord], target_anchors: List[PdbRecord]):

        """
        Constructs the PatchAligner.

        :param pfp_atoms: (list) All PdbRecord objects comprising the patch molecule.
        :param pfp_anchors: (list) PdbRecord objects of the patch atoms to be aligned.
        :param target_anchors: (list) PdbRecord objects of the target protein atoms.
        """

        self.patch_atoms = pfp_atoms  
        
        # Extract coords from provided atom objects
        pfp_anchor_coords = get_coords_array(pfp_anchors)
        target_anchor_coords = get_coords_array(target_anchors)

        pfp_centroid = np.mean(pfp_anchor_coords, axis=0)
        target_centroid = np.mean(target_anchor_coords, axis=0)

        pfp_centered = pfp_anchor_coords - pfp_centroid
        target_centered = target_anchor_coords - target_centroid

        self.rotation_object, rmsd, *_ = Rotation.align_vectors(target_centered, pfp_centered)
        rotated_pfp_centroid = self.rotation_object.apply(pfp_centroid)
        self.translation_vector = target_centroid - rotated_pfp_centroid

    def implement_align(self):
        coords = np.array([[a.x, a.y, a.z] for a in self.patch_atoms])
        rotated_coords = self.rotation_object.apply(coords)
        final_coords = rotated_coords + self.translation_vector

        for i, atom in enumerate(self.patch_atoms):
            atom.x, atom.y, atom.z = final_coords[i]
            
        return self.patch_atoms
    

    def align_single_bond(self, atoms:List[PdbRecord], at1, at2, at3, at4, target_bond_length):
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

        target_at3_pos = self.at1_coords + (target_vec_normalized * target_bond_length)

        bond_rotation_obj, rmsd = Rotation.align_vectors([-target_vec], [bond2_vec])

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

        :param p1: (tuple/list) Base parent atom (e.g., Serine CA)
        :param p2: (tuple/list) Base anchor atom (e.g., Serine N) - Axis of rotation
        :param p3: (tuple/list) Patch anchor atom (e.g., Glycan C1) - Pivot for rotation
        :param p4: (tuple/list) Patch child atom (e.g., Glycan Ring O)
        :return: (numpy.ndarray) The final transformed coordinate array.
        """
        
        self.p1_coords = np.array([p1.x, p1.y, p1.z])
        self.p2_coords = np.array([p2.x, p2.y, p2.z]) 
        self.p3_coords = np.array([p3.x, p3.y, p3.z]) 
        self.p4_coords = np.array([p4.x, p4.y, p4.z])

        current_dih = get_dihedral(self.p1_coords, self.p2_coords, self.p3_coords, self.p4_coords)
        delta = target_dih - current_dih
        rotate_dihedral(atoms, self.p3_coords, self.p2_coords, delta)
        
        return atoms

class MolGraph:

    def __init__(self, mol: Mol, distXX=1.65, distXH=1.15):

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
        self.nx_graph = None
        self.Nmols = None
        
        # Automatically build the graph upon initialization
        self._build_matrix()

    def copy(self):
        """
        Creates a shallow copy of the MolGraph. 
        """
        new_instance = copy.copy(self)
        if self.nx_graph is not None:
            new_instance.nx_graph = self.nx_graph.copy()
        return new_instance

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
            return [self.mol.records[idx].x, self.mol.records[idx].y, self.mol.records[idx].z]

        def is_hydrogen(idx):
            return self.mol.records[idx].name.strip().startswith('H')
        
        # Creates an O(1) lookup dictionary for the exact backbone Carbonyl (C) and Amide (N) atoms.
        # This allows the peptide chain to be manually linked later without global distance sweeps.
        backbone_map = {}

        spinner = ['|', '/', '-', '\\']
        for i in range(Natoms):
            # --- SPINNER (Indicates code is running)---
            if i % 30 == 0 or i == Natoms - 1:
                char = spinner[(i // 30) % 4] 
                sys.stdout.write(f'\rBuilding connectivity matrix... {char}')
                sys.stdout.flush()
            # ---------------------------------

            
            # --- MAP THE PROTEIN BACKBONE ---
            rec_i = self.mol.records[i]

            atom_name = rec_i.name.strip() 
            if atom_name in ['C', 'N']:            # Check if atom is a connecting carbonyl/amide
                key = (rec_i.chain, rec_i.res_seq) # If it is, save it to the dictionary
                if key not in backbone_map:
                    backbone_map[key] = {}
                backbone_map[key][atom_name] = i

            for j in range(i + 1, Natoms): # i + 1 prevents self-checking and double-counting
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

        sys.stdout.write('\b \n') 
        sys.stdout.flush()

        # --- CONNECT THE PROTEIN BACKBONE ---
        # Iterate through mapped residues to connect C of res n to N of res n+1
        for (chain, res_seq), atoms in backbone_map.items():
            next_key = (chain, res_seq + 1)

            if next_key in backbone_map:
                c_idx = atoms.get('C')
                n_idx = backbone_map[next_key].get('N')

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
                min_dist = float('inf')
                
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

    def __init__(self, mol_graph, moving_atoms, static_atoms, tolerance=0.75):
        self.mol_graph = mol_graph
        
        # Filter out atoms that do not move during sidechain sweeps.
        # This prevents permanent unresolvable native clashes from paralyzing the script.
        rigid_base_atoms = {'N', 'HN', 'H', 'CA', 'HA', 'C', 'O', 'CB', 'HB', 'HB1', 'HB2', 'HB3'}
        self.moving_atoms = [a for a in moving_atoms if a.name.strip() not in rigid_base_atoms]
        
        self.static_atoms = static_atoms
        self.tolerance = tolerance

        # Effective/bonded van der Waals radii (in Angstroms). See J. Charry, A. Tkatchenko, J. Chem. Theory Comput. 2024, 20, 7469–7478.
        self.vdw_radii = {
            'H': 1.10,  # modified from 1.20 since H's are artificially added in the first place
            'C': 1.70,  
            'N': 1.55,  
            'O': 1.52,  
            'P': 1.80,  
            'S': 1.80,
            'F': 1.47  # Bondi Radius. See A. Bondi, J. Phys. Chem. 1964, 68, 441–451.
        }

        self.atom_to_idx = {atom.serial: i for i, atom in enumerate(self.mol_graph.mol.records)}
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

    def check_clashes(self, limit_to_atoms=None):
        """
        Evaluates distances between moving and static atoms.
        :param limit_to_atoms: (list, optional) If provided, only check clashes originating from these PdbRecord objects.
        """
        # If no list provided, check all known moving atoms defined during init
        check_list = limit_to_atoms if limit_to_atoms is not None else self.moving_atoms

        for m_atom in check_list:
            if m_atom.name.strip() in {'N', 'H', 'CA', 'HA', 'C', 'O', 'CB'}:
                continue

            m_element = m_atom.name.strip()[0]
            m_vdw = self.vdw_radii.get(m_element, 1.50)

            for s_atom in self.static_atoms:
                if (m_atom.serial, s_atom.serial) in self.excluded_pairs:
                    continue
                if s_atom.name.strip().startswith('H'):
                    continue

                s_element = s_atom.name.strip()[0]
                s_vdw = self.vdw_radii.get(s_element, 1.50)

                actual_dist = get_distance(
                    (m_atom.x, m_atom.y, m_atom.z), 
                    (s_atom.x, s_atom.y, s_atom.z)
                )

                threshold = (m_vdw + s_vdw) * self.tolerance
                if actual_dist < threshold:
                    return f"CLASH: {m_atom.name:<4} hit {s_atom.name:<4} ({s_atom.res_seq:<3}) | Dist: {actual_dist:.2f} Å"
                    
        return None # No clash
    

class RotamerSweeper:

    """
    Generates and tests possible geometries at the patched lysine in order to resolve steric clashes.
    Sweeps through known lysine rotamers and rotates the patch around the junction where it attaches to the lysine.
    """

    def __init__(self, mol, mol_graph, res_seq, chain=" "):
        self.mol = mol
        self.graph = mol_graph
        self.res_seq = res_seq
        self.chain = chain.strip()

        self.chi_definitions = [
            ['N', 'CA', 'CB', 'CG'],   # chi1
            ['CA', 'CB', 'CG', 'CD'],  # chi2
            ['CB', 'CG', 'CD', 'CE'],  # chi3
            ['CG', 'CD', 'CE', 'NZ']   # chi4
        ]
        # Use object ID for unique mapping to avoid Serial Number duplicates
        self.obj_to_idx = {id(a): i for i, a in enumerate(self.mol.records)}

    def run_sweep(self, steric_checker):
        
        """
        Executes a three-tiered conformational search to resolve steric clashes 
        introduced by the patch alignment. 

        Tier 1 (Canonical): Tests 6 canonical staggered rotameric states for the target lysine.
        Tier 2 (Wiggle): Performs a 30-degree systematic sweep of the chi4 dihedral.
        Tier 3 (Patch Twist): Executes cumulative 15-degree twists of the patch molecule around the junction bond.

        :param steric_checker: (StericChecker) The initialized checker object to evaluate spatial overlaps.
        :return: (bool) True if a clash-free conformation is found, False if all tiers fail.
        """

        CANONICAL_LYSINE_ROTAMERS = [
            [-60, 180, 180, 180], [180, 180, 180, 180], [60, 180, 180, 180],
            [-60, -60, 180, 180], [-60, 180, 60, 180], [180, 60, 180, 180]    
        ]
        
        spinner = ['|', '/', '-', '\\']

        # --- Tier 1: Canonical Staggered Rotamers ---
        # Pivot is CA-CB. Everything from CB onwards moves.
        print("Tier 1: Trying canonical staggered rotamers...")
        t1_names = self.chi_definitions[0] # N-CA-CB-CG
        t1_pivot = self.obj_to_idx[id(self.get_atom_by_name(t1_names[1]))] # CA
        t1_axis = self.obj_to_idx[id(self.get_atom_by_name(t1_names[2]))]  # CB
        t1_moving = [self.mol.records[idx] for idx in self.get_downstream_atoms(t1_pivot, t1_axis)]
        
        for i, pose in enumerate(CANONICAL_LYSINE_ROTAMERS):
            char = spinner[i % 4]
            print(f"\rTier 1 Progress... {char}", end='', flush=True)
            
            self.apply_pose(pose) # Apply first
            clash_info = steric_checker.check_clashes(limit_to_atoms=t1_moving) # Then check
            
            if not clash_info:
                print(f"\rTier 1 Progress... Done. Success!{' ' * 50}")
                return True
            else:
                print(f"\rTier 1 Progress... {char} | {clash_info}{' ' * 10}", end='', flush=True)

        print(f"\rTier 1 Progress... Done.{' ' * 60}")
        
        # --- Tier 2: Chi4 Systematic Wiggle ---
        # Pivot is CG-CD. Only CD onwards moves.
        print("Tier 2: Canonical rotations failed. Attempting systematic chi4 sweep...")
        t2_names = self.chi_definitions[3] # CG-CD-CE-NZ
        t2_pivot = self.obj_to_idx[id(self.get_atom_by_name(t2_names[0]))] # CG
        t2_axis = self.obj_to_idx[id(self.get_atom_by_name(t2_names[1]))]  # CD
        t2_moving = [self.mol.records[idx] for idx in self.get_downstream_atoms(t2_pivot, t2_axis)]

        # Get current state of the first 3 chi angles to keep them steady
        current_pose = []
        for j in range(3):
            names = self.chi_definitions[j]
            coords = [self.get_atom_by_name(n) for n in names]
            current_pose.append(get_dihedral(*[[a.x, a.y, a.z] for a in coords]))

        for i, angle in enumerate(range(0, 360, 30)):
            char = spinner[i % 4]
            print(f"\rTier 2 Progress... {char}", end='', flush=True)
            
            self.apply_pose(current_pose + [float(angle)])
            clash_info = steric_checker.check_clashes(limit_to_atoms=t2_moving)
            
            if not clash_info:
                print(f"\rTier 2 Progress... Done. Success!{' ' * 50}")
                return True
            else:
                print(f"\rTier 2 Progress... {char} | {clash_info}{' ' * 10}", end='', flush=True)
        
        print(f"\rTier 2 Progress... Done.{' ' * 60}")

        # --- Tier 3: NZ-C7 Patch Twist ---
        print("Tier 3: Lysine stuck. Rotating patch molecule (NZ-C7)...")
        return self.attempt_patch_twist(steric_checker)

    def apply_pose(self, chi_angles):
        """Iteratively apply all chi rotations for a given pose."""
        for i, angle in enumerate(chi_angles):
            self.apply_chi_rotation(i, angle)

    def apply_chi_rotation(self, chi_index, target_angle):
        """Finds current dihedral, calculates delta, and rotates downstream atoms."""
        names = self.chi_definitions[chi_index]
        atom_records = [self.get_atom_by_name(n) for n in names]
        coords = [[a.x, a.y, a.z] for a in atom_records]
        delta = target_angle - get_dihedral(*coords)
        
        # Use object ID to look up internal graph index
        pivot_idx = self.obj_to_idx[id(atom_records[1])]
        axis_idx = self.obj_to_idx[id(atom_records[2])]

        moving_indices = self.get_downstream_atoms(pivot_idx, axis_idx)
        moving_records = [self.mol.records[i] for i in moving_indices]
        
        rotate_dihedral(moving_records, coords[1], coords[2], delta)

    def attempt_patch_twist(self, steric_checker):
        """Rotates only the patch molecule in 2-degree steps."""
        nz = self.get_atom_by_name('NZ')
        c7 = self.get_atom_by_name('C7') 
        
        pivot_idx = self.obj_to_idx[id(nz)]
        axis_idx = self.obj_to_idx[id(c7)]
        
        moving_indices_t3 = self.get_downstream_atoms(pivot_idx, axis_idx)
        moving_records_t3 = [self.mol.records[idx] for idx in moving_indices_t3]

        spinner = ['|', '/', '-', '\\']
        angles_t3 = list(range(0, 360, 2))
        
        for i, angle in enumerate(angles_t3):
            char = spinner[i % 4]
            print(f"\rTier 3 Progress... {char}", end='', flush=True)

            rotate_dihedral(moving_records_t3, (nz.x, nz.y, nz.z), (c7.x, c7.y, c7.z), 2)
            clash_info = steric_checker.check_clashes(limit_to_atoms=moving_records_t3)
            
            if not clash_info:
                print(f"\rTier 3 Progress... Done. Success!{' ' * 50}")
                return True
            else:
                print(f"\rTier 3 Progress... {char} | {clash_info}{' ' * 10}", end='', flush=True)

        print(f"\rTier 3 Progress... Done.{' ' * 60}")
        return False

    def get_atom_by_name(self, name):
        """Fetches the PdbRecord matching name, residue, AND chain."""
        try:
            # Added self.chain.strip() check here to prevent grabbing atoms from other chains
            return next(r for r in self.mol.records 
                        if r.res_seq == self.res_seq 
                        and r.name.strip() == name
                        and r.chain.strip() == self.chain)
        except StopIteration:
            print(f"Error: Could not find atom '{name}' in residue {self.res_seq} chain '{self.chain}'")
            raise

    def get_downstream_atoms(self, pivot_idx, axis_idx):
        """Uses graph connectivity to identify the side-chain segment to be moved."""
        temp_mol_graph = self.graph.copy()
        temp_nx_graph = temp_mol_graph.nx_graph

        if temp_nx_graph.has_edge(pivot_idx, axis_idx):
            temp_nx_graph.remove_edge(pivot_idx, axis_idx)
        else:
            # If the edge is missing, we can't define downstream.
            p_rec = self.mol.records[pivot_idx]
            a_rec = self.mol.records[axis_idx]
            dist = get_distance((p_rec.x, p_rec.y, p_rec.z), (a_rec.x, a_rec.y, a_rec.z))
            print(f"WARNING: Bond {pivot_idx}-{axis_idx} not found in graph. Dist: {dist:.3f}")
            return []

        # Find the component containing the axis atom
        downstream_indices = list(nx.node_connected_component(temp_nx_graph, axis_idx))
        moving_atoms = [f"{self.mol.records[i].res_name}{self.mol.records[i].res_seq}-{self.mol.records[i].name.strip()}" for i in downstream_indices]
        
        return downstream_indices