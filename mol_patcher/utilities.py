import numpy as np
import math
from scipy.spatial.distance import pdist

def get_distance(at1, at2):
    """
    Finds the distance between two atoms.
    :param at1: (list/tuple) xyz coordinates of atom1
    :param at2: (list/tuple) xyz coordinates of atom2
    :return: (float) the distance between 2 atoms
    """
    return math.sqrt((at1[0]-at2[0])**2 + (at1[1]-at2[1])**2 + (at1[2]-at2[2])**2)

def get_vector(coord1, coord2):
    """Returns the vector pointing from coord1 to coord2."""
    return np.array(coord2) - np.array(coord1)

def get_optimal_box_size(records, buffer_percent=0.33, min_buffer_nm=3.0):
    """
    Calculates the longest diagonal of the molecule and adds a dynamic box size buffer.
    The buffer scales with the molecule's size, but enforces a minimum buffer (2.0 nm) to satisfy PBC cutoffs.
    Returns the dimension and the applied buffer in nanometers.
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
    Calculates the dihedral angle (torsion angle) in degrees between four points.
    The sign follows the IUPAC convention (clockwise is positive).
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