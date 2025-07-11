import numpy as np
from numba import jit

@jit(nopython=True, cache=True, fastmath=True)
def wrap_position(position: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    """ Wrap position in a periodic lattice
    (ref: https://en.wikipedia.org/wiki/Fractional_coordinates#Relationship_between_fractional_and_Cartesian_coordinates)
    
    Args:
        position (np.ndarray): The position to wrap
        lattice (np.ndarray): The lattice of the system
    
    Returns:
        np.ndarray: The wrapped position
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
    """ Wrap positions in a periodic lattice
    (ref: https://en.wikipedia.org/wiki/Fractional_coordinates#Relationship_between_fractional_and_Cartesian_coordinates)
    
    Args:
        positions (np.ndarray): The positions to wrap
        lattice (np.ndarray): The lattice of the system
    
    Returns:
        np.ndarray: The wrapped positions
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
    """ Return the distance for a given pair of positions in a direct space."""
    return np.linalg.norm(position1 - position2)

@jit(nopython=True, cache=True, fastmath=True)
def calculate_pbc_distance(position1: np.ndarray, position2: np.ndarray, lattice: np.ndarray) -> float:
    """ Return the minimum distance between two positions, taking into account periodic boundary conditions set by the lattice.
        (ref: https://en.wikipedia.org/wiki/Fractional_coordinates#Relationship_between_fractional_and_Cartesian_coordinates)
    
    Args:
        position1 (np.ndarray): The first position
        position2 (np.ndarray): The second position
        lattice (np.ndarray): The lattice of the system
    
    Returns:
        float: The minimum distance between the two positions
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
    """ Return the angle between three positions in a direct space."""
    angle_rad = np.arccos(np.dot((position1 - position2), (position3 - position2)) / (np.linalg.norm(position1 - position2) * np.linalg.norm(position3 - position2)))
    angle_deg = np.degrees(angle_rad)
    return angle_deg

@jit(nopython=True, cache=True, fastmath=True)
def calculate_pbc_angle(position1: np.ndarray, position2: np.ndarray, position3: np.ndarray, lattice: np.ndarray) -> float:
    """ Return the angle formed by two vectors (three positions) in a periodic space.
        (ref: https://en.wikipedia.org/wiki/Fractional_coordinates#Relationship_between_fractional_and_Cartesian_coordinates)

        Args:
            position1 (np.ndarray): The first position
            position2 (np.ndarray): The second position
            position3 (np.ndarray): The third position
            lattice (np.ndarray): The lattice of the system
        
        Returns:
            float: The angle formed by the two vectors in degrees
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
    """ Convert a Cartesian position to fractional coordinates in a periodic space.
        (ref: https://en.wikipedia.org/wiki/Fractional_coordinates#Relationship_between_fractional_and_Cartesian_coordinates)
        
        Args:
            position (np.ndarray): The Cartesian position
            lattice (np.ndarray): The lattice of the system
        
        Returns:
            np.ndarray: The fractional coordinates
    """
    return np.linalg.solve(lattice.T, position.T).T

def fractional_to_cartesian(position: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    """ Convert a fractional position to Cartesian coordinates in a periodic space.
        (ref: https://en.wikipedia.org/wiki/Fractional_coordinates#Relationship_between_fractional_and_Cartesian_coordinates)
        
        Args:
            position (np.ndarray): The fractional position
            lattice (np.ndarray): The lattice of the system
        
        Returns:
            np.ndarray: The Cartesian coordinates
    """
    return np.dot(position, lattice)

__all__ = [
    'wrap_position',
    'wrap_positions',
    'calculate_direct_distance',
    'calculate_pbc_distance',
    'calculate_direct_angle',
    'calculate_pbc_angle'
]