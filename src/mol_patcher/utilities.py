import numpy as np
import math
from scipy.spatial.distance import pdist


def get_distance(at1, at2):
    """
    Calculates the Euclidean distance between two atoms in 3D space.

    :param at1: The (x, y, z) coordinates of the first atom.
    :type at1: list or tuple
    :param at2: The (x, y, z) coordinates of the second atom.
    :type at2: list or tuple
    :return: The distance between the two atoms in Angstroms.
    :rtype: float
    """

    return math.sqrt(
        (at1[0] - at2[0]) ** 2 + (at1[1] - at2[1]) ** 2 + (at1[2] - at2[2]) ** 2
    )


def get_vector(coord1, coord2):
    """
    Calculates the directional vector pointing from the first coordinate to the second.

    :param coord1: The starting (x, y, z) coordinates.
    :type coord1: list or tuple
    :param coord2: The ending (x, y, z) coordinates.
    :type coord2: list or tuple
    :return: The resulting directional vector.
    :rtype: numpy.ndarray
    """

    return np.array(coord2) - np.array(coord1)


def get_optimal_box_size(records, buffer_percent=0.33, min_buffer_nm=3.0):
    """
    Calculates the longest diagonal of the molecule and adds a dynamic box size buffer.

    The buffer scales with the molecule's size, but enforces a minimum buffer
    to satisfy periodic boundary condition (PBC) cutoffs.

    :param records: Objects representing the molecule containing x, y, z attributes.
    :type records: list
    :param buffer_percent: The percentage of the max dimension to use as a buffer, defaults to 0.33.
    :type buffer_percent: float, optional
    :param min_buffer_nm: The absolute minimum buffer size in nanometers, defaults to 3.0.
    :type min_buffer_nm: float, optional
    :return: A tuple containing the optimal box size and the applied buffer, both in nanometers.
    :rtype: tuple of (float, float)
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

    The sign follows IUPAC convention (clockwise is positive).

    :param p1: The (x, y, z) coordinates of the first atom.
    :type p1: list or tuple
    :param p2: The (x, y, z) coordinates of the second atom (axis start).
    :type p2: list or tuple
    :param p3: The (x, y, z) coordinates of the third atom (axis end).
    :type p3: list or tuple
    :param p4: The (x, y, z) coordinates of the fourth atom.
    :type p4: list or tuple
    :return: The dihedral angle in degrees (-180 to +180).
    :rtype: float
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


def identify_optimization_clusters(
    stitched_mol, base_resid, chain_id, config, base_records
):
    """
    Isolates moving versus static atoms for steric evaluation during rotamer sweeps.

    :param stitched_mol: The combined molecule object containing all records.
    :type stitched_mol: Mol
    :param base_resid: The residue sequence number of the target attachment site.
    :type base_resid: int
    :param chain_id: The chain identifier for the target residue.
    :type chain_id: str
    :param config: The configuration dictionary containing junction parameters.
    :type config: dict
    :param base_records: The original unpatched PdbRecord objects of the base protein.
    :type base_records: list
    :return: A tuple containing the list of moving atoms and the list of static atoms.
    :rtype: tuple of (list, list)
    """

    backbone_names = config.get("rigid_backbone", ["N", "CA", "C", "O"])
    is_polymer = "patch_merged_atom" not in config

    patch_res_seqs = set()
    if is_polymer:
        # Find where the patch starts in the sequence
        max_base_res = max(
            (r.res_seq for r in base_records if r.chain.strip() == chain_id.strip()),
            default=base_resid,
        )
        patch_res_seqs = set(
            r.res_seq
            for r in stitched_mol.records
            if r.chain.strip() == chain_id.strip() and r.res_seq > max_base_res
        )

    moving_atoms = [
        a
        for a in stitched_mol.records
        if (a.res_seq == base_resid and a.name.strip() not in backbone_names)
        or (a.res_seq in patch_res_seqs)
    ]

    moving_ids = {id(a) for a in moving_atoms}

    static_atoms = [a for a in stitched_mol.records if id(a) not in moving_ids]

    return moving_atoms, static_atoms
