import numpy as np
import math
from scipy.spatial.distance import pdist

def get_distance(at1, at2):
    """
    Calculates the Euclidean distance between two atoms in 3D space.

    :param at1: (list/tuple) The (x, y, z) coordinates of the first atom.
    :param at2: (list/tuple) The (x, y, z) coordinates of the second atom.
    :return: (float) The distance between the two atoms.
    """
    
    return math.sqrt((at1[0]-at2[0])**2 + (at1[1]-at2[1])**2 + (at1[2]-at2[2])**2)

def get_vector(coord1, coord2):
    """
    Calculates the directional vector pointing from the first coordinate to the second.
    
    :param coord1: (list/tuple) The starting (x, y, z) coordinates.
    :param coord2: (list/tuple) The ending (x, y, z) coordinates.
    :return: (numpy.ndarray) The resulting directional vector.
    """

    return np.array(coord2) - np.array(coord1)

def get_optimal_box_size(records, buffer_percent=0.33, min_buffer_nm=3.0):
    """
    Calculates the longest diagonal of the molecule and adds a dynamic box size buffer.
    The buffer scales with the molecule's size, but enforces a minimum buffer 
    to satisfy periodic boundary condition (PBC) cutoffs.
    
    :param records: (list) PdbRecord objects representing the molecule.
    :param buffer_percent: (float) The percentage of the max dimension to use as a buffer.
    :param min_buffer_nm: (float) The absolute minimum buffer size in nanometers.
    :return: (tuple) Optimal box size and the applied buffer, both in nanometers.
    """

    coords = np.array([[a.x, a.y, a.z] for a in records])
    
    # pdist calculates all pairwise distances
    max_dist_angstroms = np.max(pdist(coords))
    max_dist_nm = max_dist_angstroms / 10.0
    
    # Calculate the percentage-based buffer
    calculated_buffer = max_dist_nm * buffer_percent
    
    # Enforce the minimum
    final_buffer = max(calculated_buffer, min_buffer_nm)
    
    # Total size and round to 3 decimal places
    optimal_size = round(max_dist_nm + final_buffer, 3)
    
    return optimal_size, round(final_buffer, 3)

def get_dihedral(p1, p2, p3, p4):
    """
    Calculates the dihedral angle in degrees between four points.
    The sign follows the IUPAC convention (clockwise is positive).
    
    :param p1: (list/tuple) The (x, y, z) coordinates of the first atom.
    :param p2: (list/tuple) The (x, y, z) coordinates of the second atom (axis start).
    :param p3: (list/tuple) The (x, y, z) coordinates of the third atom (axis end).
    :param p4: (list/tuple) The (x, y, z) coordinates of the fourth atom.
    :return: (float) The dihedral angle in degrees (-180 to +180).
    """
    # Create vectors between the atoms
    b0 = -1.0 * (np.array(p2) - np.array(p1))
    b1 = np.array(p3) - np.array(p2)
    b2 = np.array(p4) - np.array(p3)

    # Normalize the central bond (the axis of rotation)
    b1 /= np.linalg.norm(b1)

    # Find the normal vectors to the two planes
    # v = vector in the first plane (A-B-C)
    # w = vector in the second plane (B-C-D)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    # Calculate the angle using the dot product and the determinant
    # x is the cosine component, y is the sine component
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)

    # Return the angle in degrees
    return np.degrees(np.arctan2(y, x))

def identify_optimization_clusters(stitched_mol, base_resid, chain_id, config, base_records):
    """Isolates moving vs static atoms for the StericChecker."""
    backbone_names = config.get("rigid_backbone", ['N', 'CA', 'C', 'O'])
    is_polymer = "patch_merged_atom" not in config
    
    patch_res_seqs = set()
    if is_polymer:
        # Find where the patch starts in the sequence
        max_base_res = max((r.res_seq for r in base_records if r.chain.strip() == chain_id.strip()), default=base_resid)
        patch_res_seqs = set(r.res_seq for r in stitched_mol.records 
                            if r.chain.strip() == chain_id.strip() and r.res_seq > max_base_res)

    moving_atoms = [a for a in stitched_mol.records 
                    if (a.res_seq == base_resid and a.name.strip() not in backbone_names) 
                    or (a.res_seq in patch_res_seqs)]
    static_atoms = [a for a in stitched_mol.records if a not in moving_atoms]
    
    return moving_atoms, static_atoms