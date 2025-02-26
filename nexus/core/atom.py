"""
Atom module for the Nexus project.

This module defines the Atom class, which represents an atom within a system.

Imports:
    - Standard libraries
    - Third-party libraries
    - Internal modules

Classes:
    - Atom: Represents an atom within a system.
"""

# Standard library imports
import sys

# Third-party imports
import numpy as np

# Internal imports
from ..data import chemical_symbols, atomic_masses
from .cutoff import Cutoff

class Atom:
    """
    Represents an atom within a system.

    Attributes:
        element (str): Atomic element.
        id (int): Identifier of the atom in the system.
        position (np.array): XYZ coordinates.
        frame (int): Frame index in the trajectory.
        cutoffs (dict): Cutoff distances dictionary (Cutoff object).
        extension (str): Extension used for method determination.
        atomic_mass (float): Atomic mass of the atom.
        neighbours (list): List of first neighbours (PBC applied).
        coordination (int): Number of neighbours around the atom (PBC applied).
        parent (Atom): Root of the cluster.
        cluster_id (int): Cluster index to which the atom belongs.
        cluster_type (str): Cluster type to which the atom belongs.

    Methods:
        __init__: Initializes an Atom object.
        get_element: Returns the element of the Atom.
        get_id: Returns the unique identifier of the Atom.
        get_position: Returns the spatial coordinates of the Atom.
        get_frame: Returns the frame index associated with the Atom.
        get_neighbours: Returns the list of neighbours of the Atom.
        get_atomic_mass: Returns the atomic mass of the Atom.
        get_coordination: Returns the coordination number of the Atom.
        get_cluster_id: Returns the cluster id that the atom belongs to.
        add_neighbour: Adds a neighbour to the list of neighbours of the Atom.
        filter_neighbours: Removes neighbours not within cutoff distances (depending on pair of atoms).
        set_cluster: Sets the cluster and type that the atom belongs to.
        reset_cluster: Resets the cluster(s) of the atom.
    """
    def __init__(self, element, id, position, frame, cutoffs, extension="SiOz") -> None:
        """
        Initializes an Atom object with the provided information.

        Parameters:
            element (str): Atomic element.
            id (int): Identifier of the atom in the system.
            position (np.array): XYZ coordinates.
            frame (int): Frame index in the trajectory.
            cutoffs (dict): Cutoff distances dictionary (Cutoff object).
            extension (str): Extension used for method determination.
        """
        # Initialize an Atom object with the provided information
        self.element : str = element                    # atomic element
        self.id : int = id                              # id of the atom in the system
        self.position : np.array = np.array(position)   # xyz coordinates
        self.frame : int = frame                        # frame that this atom belong to in the trajectory
        self.cutoffs : Cutoff = cutoffs                 # cutoffs dictionary (Cutoff object) 
        
        # Initialize the extension so that correct methods are used.
        self.extension : str = extension
        
        # Initialize atomic data from the periodic table and other informations
        if self.element in chemical_symbols:
            index = np.where(self.element == chemical_symbols)[0].astype(int)
            self.atomic_mass : float = atomic_masses[index][0]
        else:
            print(f"\tERROR: Element {self.element} not found in the periodic table.")
            print("\tFailed to initialize the Atom object {self.id} in the frame {self.frame}.")
            print("Exiting.")
            sys.exit(1)
        
        # Initialize neighbours attributes 
        self.neighbours : list = [] # first neighbours (pbc applied)
        self.coordination : int = 0 # number of neighbours around the atom (pbc applied)
        
        # Initialize some attributes for clustering 
        self.parent : Atom = self       # parent is the root of the cluster. 
        self.cluster_id : int = None    # cluster index in which the atom belong to
        self.cluster_type : str = None  # cluster type in which the atom belong to

    #____________GETTER METHODS____________
    
    def get_element(self) -> str:
        """
        Return the element of the Atom.
        
        Returns:
            str : Name of the element of the Atom.
        """
        return self.element
    
    def get_id(self) -> int:
        """
        Return the unique identifier of the Atom.
        
        Returns:
            int : Index of the Atom in the frame.
        """
        return self.id
    
    def get_position(self) -> np.array:
        """
        Return the spatial coordinates of the Atom.
        
        Returns:
            np.array : Cartesian coordinates of the Atom.
        """
        return self.position
    
    def get_frame(self) -> int:
        """
        Return the frame index associated with the Atom.
        
        Returns:
            int : Frame index of the trajectory.
        """
        return self.frame
    
    def get_neighbours(self) -> list:
        """
        Return the list of neighbours of the Atom.
        
        Returns:
            list : List of the nearest neighbours of the Atom.
        """
        return self.neighbours
    
    def get_atomic_mass(self) -> float:
        """
        Return the atomic mass of the Atom.
        
        Returns:
            float : Atomic mass of the Atom.
        """
        return self.atomic_mass

    def get_coordination(self) -> int:
        """
        Return the coordination number of the Atom. (ie the number of first neighbours)
        
        Returns:
            int : Coordination number of the Atom.
        """
        return self.coordination

    def get_cluster_id(self) -> int:
        """
        Return the cluster id that the atom belong to.
        
        Returns:
            int : Cluster index that the Atom belong to.
        """
        return self.cluster_id
    
    #____________NEIGHBOURS METHODS____________
    
    def add_neighbour(self, neighbour) -> None:
        """
        Add a neighbour to the list of neighbours of the Atom.
        
        Parameters:
            neighbour (Atom) : Atom object to append to the list of neighbours.
            
        Returns:
            None.
        """
        self.neighbours.append(neighbour)
    
    def filter_neighbours(self, distances) -> None:
        """
        Removes the neighbours that are not within the cutoff distances (depending of pair of atoms).
        
        Parameters:
            distances (list): List of distances to the neighbours.
        
        Returns:
            None.
        """
        new_list_neighbours = []
        new_list_distances  = []
        
        for k, neighbour in enumerate(self.neighbours):
            rcut = self.cutoffs.get_cutoff(self.element, neighbour.get_element())
            
            if isinstance(distances, float):
                # if 'distances' is a float, it means that the neighbour of this atom is itself.
                current_distance = distances
            else:
                current_distance = distances[k]
            
            if current_distance > rcut: # neighbour is too far 
                continue # go next neighbour
            elif current_distance == 0: # neighbour is this atom.
                continue # go next neighbour
            else:
                new_list_neighbours.append(neighbour) # keep the neighbour
                new_list_distances.append(current_distance)

        self.neighbours = new_list_neighbours
    
    #____________CLUSTERS METHODS____________
    
    def set_cluster(self, cluster_id, cluster_type) -> None:
        """
        Set the cluster if and type that the atom belong to.
        
        Parameters:
            cluster_id (int) : cluster index to be set.
            cluster_type (str) : cluster type to be set. 
        
        Returns:
            None.
        """
        self.cluster_id = cluster_id
        self.cluster_type = cluster_type
    
    def reset_cluster(self) -> None:
        """
        Reset the cluster(s) of the atom.
        
        Returns:
            None.
        """
        self.cluster_id     = None
        self.cluster_type   = None
        self.parent         = self
