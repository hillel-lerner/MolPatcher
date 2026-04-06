import numpy as np
from scipy.spatial.transform import Rotation
from typing import List, Tuple
from .mol_record import PdbRecord, Mol
from .utilities import get_distance, get_dihedral
from .mol_record import Mol
import networkx as nx
import math
import sys

    
def get_coords_array(atoms: List[PdbRecord]) -> np.ndarray:
    """
    Extracts an Nx3 numpy array of coordinates from any list of PdbRecord objects.
    """
    return np.array([[a.x, a.y, a.z] for a in atoms])

def rotate_dihedral(atoms, pivot_coord, axis_coord, angle_degrees):
        """
        Applies a quaternion rotation to a list of atoms around a specific bond vector.
        
        :param atoms: List of PdbRecord objects to be rotated.
        :param pivot_coord: (x,y,z) of the stationary atom (e.g., CA).
        :param axis_coord: (x,y,z) of the atom defining the rotation axis (e.g., CB).
        :param angle_degrees: The amount to rotate in degrees.
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
        for atom in atoms:
            pos = np.array([atom.x, atom.y, atom.z])
            
            # Translate atom so pivot is at origin, rotate, then translate back
            new_pos = rot_obj.apply(pos - pivot) + pivot

            # Update the PdbRecord attributes
            atom.x, atom.y, atom.z = new_pos
        return new_pos

class PatchAligner:
    """
    Takes a list of EQUIVALENT PdbRecord objects from a patch molecule (pfp) and base molecule (protein) as anchors
    Generates the rotational and transformational matrices required to align the anchors
    Apply those operations to all atoms of patch molecule
    Update the coordinates of the patch molecules PdbRecord Objects
    """

    def __init__(self, pfp_atoms: List[PdbRecord], pfp_anchors: List[PdbRecord], target_anchors: List[PdbRecord]):
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

    def _build_matrix(self):
        """Internal method: executes distance checks and builds the connectivity matrix."""
        Natoms = len(self.mol.records)
        conn_mat = np.zeros((Natoms, Natoms))


        def get_coords(idx):
            return [self.mol.records[idx].x, self.mol.records[idx].y, self.mol.records[idx].z]

        def is_hydrogen(idx):
            return self.mol.records[idx].name.strip().startswith('H')

        spinner = ['|', '/', '-', '\\']
        for i in range(Natoms):
            # --- SPINNER (Indicates code is running)---
            if i % 30 == 0 or i == Natoms - 1:
                char = spinner[(i // 30) % 4] 
                sys.stdout.write(f'\rBuilding connectivity matrix... {char}')
                sys.stdout.flush()
            # ---------------------------------

            for j in range(i + 1, Natoms): # i + 1 prevents self-checking and double-counting
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
        self.moving_atoms = moving_atoms
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

    def check_clashes(self):
        for m_atom in self.moving_atoms:
            
            m_element = m_atom.name.strip()[0]
            m_vdw = self.vdw_radii.get(m_element, 1.50)

            for s_atom in self.static_atoms:
                
                # Check exclusion list FIRST (O(1) lookup)
                if (m_atom.serial, s_atom.serial) in self.excluded_pairs:
                    continue

                s_element = s_atom.name.strip()[0]
                s_vdw = self.vdw_radii.get(s_element, 1.50)

                # Calculate the 3D distance
                actual_dist = get_distance(
                    (m_atom.x, m_atom.y, m_atom.z), 
                    (s_atom.x, s_atom.y, s_atom.z)
                )

                # Check for steric clash
                threshold = (m_vdw + s_vdw) * self.tolerance
                if actual_dist < threshold:
                    print(f"CLASH: {m_atom.name} ({m_atom.res_seq}) hit {s_atom.name} ({s_atom.res_seq}) | Dist: {actual_dist:.2f}")
                    return True 
                    
        return False
    

class RotamerSweeper:

    """
    Coordinates side-chain and patch rotations to find a valid protein conformation.
    Uses a three-tiered search strategy (Golden, Wiggle, and Patch-Tweak).
    """

    def __init__(self, mol, mol_graph, res_seq):
        """
        Initializes the sweeper for a specific lysine.
        
        :param mol: The Mol object containing records and atoms.
        :param mol_graph: The pre-built MolGraph (connectivity matrix).
        :param res_seq: The residue number to be rotated (i.e., target Lysine).
        """
        self.mol = mol
        self.graph = mol_graph
        self.res_seq = res_seq
        self.chi_definitions = [
            ['N', 'CA', 'CB', 'CG'],   # chi1
            ['CA', 'CB', 'CG', 'CD'],  # chi2
            ['CB', 'CG', 'CD', 'CE'],  # chi3
            ['CG', 'CD', 'CE', 'NZ']   # chi4
        ]
        self.serial_to_idx = {a.serial: i for i, a in enumerate(self.mol.records)}


    def run_sweep(self, steric_checker):
        """Executes the three-tiered conformational search strategy."""
        GOLDEN_LYSINE_ROTAMERS = [
            [-60, 180, 180, 180], [180, 180, 180, 180], [60, 180, 180, 180],
            [-60, -60, 180, 180], [-60, 180, 60, 180], [180, 60, 180, 180]    
        ]
        
        spinner = ['|', '/', '-', '\\']

        # Tier 1: Canonical Staggered Rotamers 
        print("Tier 1: Trying canonical staggered rotamers...")
        for i, pose in enumerate(GOLDEN_LYSINE_ROTAMERS):
            sys.stdout.write(f'\rTier 1 Progress... {spinner[i % 4]}')
            sys.stdout.flush()
            
            self.apply_pose(pose)
            if not steric_checker.check_clashes(): 
                sys.stdout.write('\b Success!\n')
                return True
        sys.stdout.write('\b Done.\n')
        
        # Tier 2: Chi4 Systematic Wiggle 
        print("Tier 2: Golden failed. Attempting systematic chi4 sweep...")
        base_pose = [-60, 180, 180] 
        angles_t2 = list(range(0, 360, 30))
        
        for i, angle in enumerate(angles_t2):
            sys.stdout.write(f'\rTier 2 Progress... {spinner[i % 4]}')
            sys.stdout.flush()
            
            self.apply_pose(base_pose + [angle])
            if not steric_checker.check_clashes(): 
                sys.stdout.write('\b Success!\n')
                return True
        sys.stdout.write('\b Done.\n')

        # Tier 3: NZ-C7 Patch Tweak
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
        
        # Pivot (CA) stays still, Axis (CB) defines the rotation vector
        moving_atoms = self.get_downstream_atoms(self.serial_to_idx[atom_records[1].serial], 
                                            self.serial_to_idx[atom_records[2].serial])
        rotate_dihedral(moving_atoms, coords[1], coords[2], delta)

    def attempt_patch_twist(self, steric_checker):
        """Rotates only the patch molecule in 15-degree steps."""
        nz = self.get_atom_by_name('NZ')
        c7 = self.get_atom_by_name('C7') 
        moving_atoms = self.get_downstream_atoms(self.serial_to_idx[nz.serial], self.serial_to_idx[c7.serial])

        spinner = ['|', '/', '-', '\\']
        angles_t3 = list(range(0, 360, 15))
        
        for i, angle in enumerate(angles_t3):
            sys.stdout.write(f'\rTier 3 Progress... {spinner[i % 4]}')
            sys.stdout.flush()
            
            rotate_dihedral(moving_atoms, (nz.x, nz.y, nz.z), (c7.x, c7.y, c7.z), 15)
            if not steric_checker.check_clashes():
                sys.stdout.write('\b Success!\n')
                print(f"Success! Patch rotated {angle+15} degrees to clear clash.")
                return True
                
        sys.stdout.write('\b Done.\n')
        return False

    def get_atom_by_name(self, name):
        """Fetches the PdbRecord matching the name within the target residue."""
        return next(r for r in self.mol.records if r.res_seq == self.res_seq and r.name.strip() == name)

    def get_downstream_atoms(self, pivot_idx, axis_idx):
        """Uses graph connectivity to identify the side-chain segment to be moved."""
        temp_graph = self.graph.nx_graph.copy()
        temp_graph.remove_edge(pivot_idx, axis_idx)
        moving_indices = nx.node_connected_component(temp_graph, axis_idx)
        return [self.mol.records[i] for i in moving_indices]