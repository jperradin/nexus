import numpy as np
from numba import jit

@jit(nopython=True, cache=True, fastmath=True)
def wrap_position(position: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    """
    Wrap a single position into the simulation cell using fractional coordinates.

    Args:
        position (np.ndarray): Cartesian position vector of length 3.
        lattice (np.ndarray): 3x3 lattice matrix defining the simulation cell.

    Returns:
        np.ndarray: Wrapped Cartesian position inside the cell.
    """
    position = np.ascontiguousarray(position)
    lattice = np.ascontiguousarray(lattice)

    lattice_inv = np.linalg.inv(lattice)
    lattice_inv = np.ascontiguousarray(lattice_inv)

    fractional_position = np.dot(position, lattice_inv)
    fractional_position = np.ascontiguousarray(fractional_position)

    fractional_position -= np.floor(fractional_position)
    fractional_position = np.ascontiguousarray(fractional_position)

    wrapped_position = np.dot(fractional_position, lattice)
    wrapped_position = np.ascontiguousarray(wrapped_position)

    return wrapped_position

@jit(nopython=True, cache=True, fastmath=True)
def wrap_positions(positions: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    """
    Wrap an array of positions into the simulation cell using fractional coordinates.

    Args:
        positions (np.ndarray): Cartesian positions of shape (N, 3).
        lattice (np.ndarray): 3x3 lattice matrix defining the simulation cell.

    Returns:
        np.ndarray: Wrapped Cartesian positions of shape (N, 3).
    """
    wrapped_positions = np.zeros_like(positions)
    for i in range(positions.shape[0]):
        position = np.ascontiguousarray(positions[i])
        lattice = np.ascontiguousarray(lattice)

        lattice_inv = np.linalg.inv(lattice)
        lattice_inv = np.ascontiguousarray(lattice_inv)

        fractional_position = np.dot(position, lattice_inv)
        fractional_position = np.ascontiguousarray(fractional_position)

        fractional_position -= np.floor(fractional_position)
        fractional_position = np.ascontiguousarray(fractional_position)

        wrapped_position = np.dot(fractional_position, lattice)
        wrapped_position = np.ascontiguousarray(wrapped_position)

        wrapped_positions[i] = wrapped_position

    return wrapped_positions

@jit(nopython=True, cache=True, fastmath=True)
def calculate_direct_distance(position1: np.ndarray, position2: np.ndarray) -> float:
    """
    Compute the Euclidean distance between two positions in direct space.

    Args:
        position1 (np.ndarray): First Cartesian position vector.
        position2 (np.ndarray): Second Cartesian position vector.

    Returns:
        float: Euclidean distance between the two positions.
    """
    return np.linalg.norm(position1 - position2)

@jit(nopython=True, cache=True, fastmath=True)
def calculate_pbc_distance(position1: np.ndarray, position2: np.ndarray, lattice: np.ndarray) -> float:
    """
    Compute the minimum-image distance between two positions under periodic boundary
    conditions.

    Converts the displacement to fractional coordinates, applies the minimum-image
    convention via rounding, and converts back to Cartesian space.

    Args:
        position1 (np.ndarray): First Cartesian position vector.
        position2 (np.ndarray): Second Cartesian position vector.
        lattice (np.ndarray): 3x3 lattice matrix defining the simulation cell.

    Returns:
        float: Minimum-image distance between the two positions.
    """
    # Calculate the direct displacement vector
    direct_disp = position1 - position2

    # Get fractional coordinates in the lattice 
    inv_lattice = np.linalg.inv(lattice)
    frac_disp = np.dot(inv_lattice, direct_disp)

    # Minimum image convention
    frac_disp -= np.round(frac_disp)

    min_disp = np.dot(frac_disp, lattice)

    return np.linalg.norm(min_disp)

@jit(nopython=True, cache=True, fastmath=True)
def calculate_direct_angle(position1: np.ndarray, position2: np.ndarray, position3: np.ndarray) -> float:
    """
    Compute the angle formed by three positions in direct space.

    The angle is measured at ``position2`` between the vectors to ``position1`` and
    ``position3``.

    Args:
        position1 (np.ndarray): First endpoint Cartesian position.
        position2 (np.ndarray): Vertex Cartesian position.
        position3 (np.ndarray): Second endpoint Cartesian position.

    Returns:
        float: Angle in degrees.
    """
    angle_rad = np.arccos(np.dot((position1 - position2), (position3 - position2)) / (np.linalg.norm(position1 - position2) * np.linalg.norm(position3 - position2)))
    angle_deg = np.degrees(angle_rad)
    return angle_deg

@jit(nopython=True, cache=True, fastmath=True)
def calculate_pbc_angle(position1: np.ndarray, position2: np.ndarray, position3: np.ndarray, lattice: np.ndarray) -> float:
    """
    Compute the angle formed by three positions under periodic boundary conditions.

    The angle is measured at ``position2``. Displacement vectors are converted to
    fractional coordinates and wrapped via the minimum-image convention before computing
    the angle.

    Args:
        position1 (np.ndarray): First endpoint Cartesian position.
        position2 (np.ndarray): Vertex Cartesian position.
        position3 (np.ndarray): Second endpoint Cartesian position.
        lattice (np.ndarray): 3x3 lattice matrix defining the simulation cell.

    Returns:
        float: Angle in degrees.
    """
    # Calculate the direct displacement vectors assuming position2 is the center
    direct_disp1 = position1 - position2
    direct_disp2 = position2 - position3

    # Get fractional coordinates in the lattice 
    inv_lattice = np.linalg.inv(lattice)
    frac_disp1 = np.dot(inv_lattice, direct_disp1)
    frac_disp2 = np.dot(inv_lattice, direct_disp2)

    # Minimum image convention
    frac_disp1 -= np.round(frac_disp1)
    frac_disp2 -= np.round(frac_disp2)

    min_disp1 = np.dot(frac_disp1, lattice)
    min_disp2 = np.dot(frac_disp2, lattice)

    # Calculate the angle
    dot_product = np.dot(min_disp1, min_disp2)
    norm1 = np.linalg.norm(min_disp1)
    norm2 = np.linalg.norm(min_disp2)

    # cos_angle = np.clip(dot_product / (norm1 * norm2), -1.0, 1.0)
    cos_angle = dot_product / (norm1 * norm2)
    cos_angle = max(-1.0, min(1.0, cos_angle)) # equivalent to np.clip

    angle_rad = np.arccos(cos_angle)
    angle_deg = np.degrees(angle_rad)

    return angle_deg

def cartesian_to_fractional(position: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    """
    Convert Cartesian coordinates to fractional coordinates.

    Args:
        position (np.ndarray): Cartesian position(s) to convert.
        lattice (np.ndarray): 3x3 lattice matrix defining the simulation cell.

    Returns:
        np.ndarray: Fractional coordinates in the lattice basis.
    """
    return np.linalg.solve(lattice.T, position.T).T

def fractional_to_cartesian(position: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    """
    Convert fractional coordinates to Cartesian coordinates.

    Args:
        position (np.ndarray): Fractional position(s) to convert.
        lattice (np.ndarray): 3x3 lattice matrix defining the simulation cell.

    Returns:
        np.ndarray: Cartesian coordinates.
    """
    return np.dot(position, lattice)

@jit(nopython=True, cache=True, fastmath=True)
def calculate_gyration_radius(positions: np.ndarray, center_of_mass: np.ndarray) -> float:
    """
    Compute the radius of gyration for a set of positions around a center of mass.

    Calculates the root-mean-square distance of all positions from the center of mass.

    Args:
        positions (np.ndarray): Cartesian positions of shape (N, 3).
        center_of_mass (np.ndarray): Center-of-mass position vector of length 3.

    Returns:
        float: Radius of gyration, or 0.0 if positions is empty.
    """
    if positions.shape[0] == 0:
        return 0.0

    rg_squared = 0.0
    n_nodes = positions.shape[0]

    # Sum the squared distances from the center of mass
    for i in range(n_nodes):
        dx = positions[i, 0] - center_of_mass[0]
        dy = positions[i, 1] - center_of_mass[1]
        dz = positions[i, 2] - center_of_mass[2]
        rg_squared += dx**2 + dy**2 + dz**2

    # Return the root of the mean squared distance
    return np.sqrt(rg_squared / n_nodes)    

@jit(nopython=True, cache=True, fastmath=True)
def filter_neighbors_pbc_batch(
    target: np.ndarray,
    cands: np.ndarray,
    lattice: np.ndarray,
    inv_lattice: np.ndarray,
    rcut2: np.ndarray,
):
    """
    Vectorized PBC-distance filter for one node against many candidates.

    Args:
        target (np.ndarray): (3,) cartesian position of the central node.
        cands  (np.ndarray): (M, 3) cartesian positions of candidates.
        lattice (np.ndarray): (3, 3) lattice matrix.
        inv_lattice (np.ndarray): (3, 3) precomputed inverse lattice.
        rcut2 (np.ndarray): (M,) per-pair squared cutoff. Negative entries flag
            invalid pairs and are skipped.

    Returns:
        keep (np.ndarray): (M,) bool mask of accepted candidates.
        dists (np.ndarray): (M,) float64, distance for kept entries (0 elsewhere).
    """
    M = cands.shape[0]
    keep = np.zeros(M, dtype=np.bool_)
    dists = np.zeros(M, dtype=np.float64)
    tx = target[0]; ty = target[1]; tz = target[2]
    for i in range(M):
        rc2 = rcut2[i]
        if rc2 < 0.0:
            continue
        dx = cands[i, 0] - tx
        dy = cands[i, 1] - ty
        dz = cands[i, 2] - tz
        # Row convention: frac = d @ inv_lattice, cart = frac @ lattice
        fx = inv_lattice[0, 0] * dx + inv_lattice[1, 0] * dy + inv_lattice[2, 0] * dz
        fy = inv_lattice[0, 1] * dx + inv_lattice[1, 1] * dy + inv_lattice[2, 1] * dz
        fz = inv_lattice[0, 2] * dx + inv_lattice[1, 2] * dy + inv_lattice[2, 2] * dz
        fx -= np.round(fx); fy -= np.round(fy); fz -= np.round(fz)
        mx = lattice[0, 0] * fx + lattice[1, 0] * fy + lattice[2, 0] * fz
        my = lattice[0, 1] * fx + lattice[1, 1] * fy + lattice[2, 1] * fz
        mz = lattice[0, 2] * fx + lattice[1, 2] * fy + lattice[2, 2] * fz
        d2 = mx * mx + my * my + mz * mz
        if d2 <= rc2:
            keep[i] = True
            dists[i] = np.sqrt(d2)
    return keep, dists


@jit(nopython=True, cache=True, fastmath=True)
def filter_neighbors_direct_batch(
    target: np.ndarray,
    cands: np.ndarray,
    rcut2: np.ndarray,
):
    """
    Vectorized direct-distance (no PBC) filter, otherwise same contract as
    ``filter_neighbors_pbc_batch``.

    Args:
        target (np.ndarray): (3,) cartesian position of the central node.
        cands  (np.ndarray): (M, 3) cartesian positions of candidates.
        rcut2 (np.ndarray): (M,) per-pair squared cutoff. Negative entries flag
            invalid pairs and are skipped.

    Returns:
        keep (np.ndarray): (M,) bool mask of accepted candidates.
        dists (np.ndarray): (M,) float64, distance for kept entries (0 elsewhere).
    """
    M = cands.shape[0]
    keep = np.zeros(M, dtype=np.bool_)
    dists = np.zeros(M, dtype=np.float64)
    tx = target[0]; ty = target[1]; tz = target[2]
    for i in range(M):
        rc2 = rcut2[i]
        if rc2 < 0.0:
            continue
        dx = cands[i, 0] - tx
        dy = cands[i, 1] - ty
        dz = cands[i, 2] - tz
        d2 = dx * dx + dy * dy + dz * dz
        if d2 <= rc2:
            keep[i] = True
            dists[i] = np.sqrt(d2)
    return keep, dists


__all__ = [
    'wrap_position',
    'wrap_positions',
    'calculate_direct_distance',
    'calculate_pbc_distance',
    'calculate_direct_angle',
    'calculate_pbc_angle',
    'calculate_gyration_radius',
    'filter_neighbors_pbc_batch',
    'filter_neighbors_direct_batch',
]
