import numpy as np
from scipy.spatial.transform import Rotation
from typing import List, Tuple
from .mol_record import PdbRecord

class FindAnchors:
    @staticmethod
    def extract_coords(pfp_anchors: List[PdbRecord], target_anchors: List[PdbRecord]) -> Tuple[np.ndarray, np.ndarray]:
        # pfp_anchors are hardcoded in main.py
        # target anchor atoms are listed in main.py 
        pfp_coords = [[a.x, a.y, a.z] for a in pfp_anchors]
        target_coords = [[a.x, a.y, a.z] for a in target_anchors]
        return np.array(pfp_coords), np.array(target_coords)
    

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
        pfp_anchor_coords, target_anchor_coords = FindAnchors.extract_coords(pfp_anchors, target_anchors)

        pfp_centroid = np.mean(pfp_anchor_coords, axis=0)
        target_centroid = np.mean(target_anchor_coords, axis=0)

        pfp_centered = pfp_anchor_coords - pfp_centroid
        target_centered = target_anchor_coords - target_centroid

        self.rotation_object, rmsd = Rotation.align_vectors(target_centered, pfp_centered)
        rotated_pfp_centroid = self.rotation_object.apply(pfp_centroid)
        self.translation_vector = target_centroid - rotated_pfp_centroid

    def implement_align(self):
        coords = np.array([[a.x, a.y, a.z] for a in self.patch_atoms])
        rotated_coords = self.rotation_object.apply(coords)
        final_coords = rotated_coords + self.translation_vector

        for i, atom in enumerate(self.patch_atoms):
            atom.x, atom.y, atom.z = final_coords[i]
            
        return self.patch_atoms