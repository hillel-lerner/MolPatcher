import numpy as np
from scipy.spatial.transform import Rotation
from typing import List, Tuple
from .mol_record import PdbRecord, Mol
from .utilities import get_distance
import networkx as nx
import math
from .mol_record import Mol
    
def get_coords_array(atoms: List[PdbRecord]) -> np.ndarray:
    """
    Extracts an Nx3 numpy array of coordinates from any list of PdbRecord objects.
    """
    return np.array([[a.x, a.y, a.z] for a in atoms])

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



        for i in range(Natoms):
            for j in range(i + 1, Natoms): # i + 1 prevents self-checking and double-counting
                dist = get_distance(get_coords(i), get_coords(j))
                
                is_h = is_hydrogen(i) or is_hydrogen(j)

                if is_h and dist <= self.distXH:
                    conn_mat[i, j] = 1
                    conn_mat[j, i] = 1

                elif not is_h and dist <= self.distXX:
                    conn_mat[i, j] = 1
                    conn_mat[j, i] = 1  

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