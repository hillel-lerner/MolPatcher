import numpy as np
import math

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