"""
Cluster module for the Nexus project.

This module defines the Cluster class, which represents a cluster of atoms within a system.

Imports:
    - Standard libraries
    - Third-party libraries
    - Internal modules

Classes:
    - Cluster: Represents a cluster of atoms within a system.
"""

# Standard library imports
import os

# Third-party imports
import numpy as np
from tqdm import tqdm

# Internal imports
from .box import Box

class Cluster:
    """
    Represents a cluster of atoms within a system.

    Attributes:
        atoms (list): List of Atom objects belonging to the cluster.
        box (Box): Box object representing the simulation box.
        connectivity (str): Connectivity of the cluster.
        root_id (int): Atom id that is the root of the cluster.
        frame (int): Frame of the trajectory.
        size (int): Size of the cluster (number of atoms).
        center_of_mass (list): Center of mass of the cluster.
        indices (list): List of indices of the atoms.
        unwrapped_positions (list): List of unwrapped positions of the cluster.
        percolation_probability (str): Percolation probability.
        gyration_radius (float): Gyration radius of the cluster.

    Methods:
        __init__: Initializes a Cluster object.
        add_atom: Adds an atom to the cluster.
        get_atoms: Returns the list of Atom objects belonging to the cluster.
        set_indices_and_positions: Sets the array of unique indices and positions of atoms in the cluster.
        calculate_center_of_mass: Calculates the center of mass of the cluster.
        write_coordinates: Writes the cluster coordinates to an XYZ file.
        calculate_gyration_radius: Calculates the gyration radius of the cluster.
        calculate_percolation_probability: Calculates the percolation probability of the cluster.
        calculate_unwrapped_positions: Calculates the unwrapped positions of atoms in the cluster.
        unwrap_position: Unwraps position considering periodic boundary conditions.
    """
    def __init__(self, atoms=None, box=None, connectivity="", root_id=None, frame=None, size=None) -> None:
        """
        Initializes a Cluster object.

        Parameters:
            atoms (list): List of Atom objects belonging to the cluster.
            box (Box): Box object representing the simulation box.
            connectivity (str): Connectivity of the cluster.
            root_id (int): Atom id that is the root of the cluster.
            frame (int): Frame of the trajectory.
            size (int): Size of the cluster (number of atoms).
        """
        # Initialize Cluster object with the provided information in argument.
        self.atoms : list = atoms if atoms is not None else []  # list of Atom objects belonging to the cluster
        self.box : Box = box                                    # Box object
        self.connectivity : str = connectivity                  # type of connectivity of the cluster
        self.root_id : int = root_id                            # Atom id that is the root of cluster
        self.frame : int = frame                                # frame of the trajectory
        self.size : int = size                                  # size of cluster (number of atoms)
        
        # Attributes that has to be calculated.
        self.center_of_mass : list = []         # Center of mass of the cluster
        self.indices : list = []                # List of indices of the Atoms
        self.unwrapped_positions : list = []    # List of unwrapped positions of the cluster
        self.percolation_probability : str = '' # Percolation probability, ie 'x' the percolates through the x dimension
        self.gyration_radius : float = 0.0      # Gyration radius of the cluster
        self.order_parameter : list = [0.0] * 3 # Order parameter of the cluster [1d value, 2d value, 3d value]
        self.number_of_nodes : int = 0          # Number of nodes in all clusters of the same connectivity
        self.concentration : float = 0.0        # Concentration of the 
        
    def add_atom(self, atom) -> None:
        """
        Add an atom to the cluster.
        
        Parameters:
            atom (Atom): Atom object to be added to the cluster
        """
        self.atoms.append(atom)
        
    def get_atoms(self) -> list:
        """
        Return the list of Atom objects belonging to the cluster
        
        Returns:
            None.
        """
        return self.atoms
    
    def set_indices_and_positions(self, positions_dict) -> None:
        """
        Set the array of all unique indices of Atom objects in the cluster.

        Parameters:
            positions_dict (dict): Dictionary where keys are the index of the Atoms and values, their unwrapped positions.
        
        Returns:
            None.
        """
        for atom_id, position in positions_dict.items():
            self.indices.append(atom_id)
            self.unwrapped_positions.append([position[0], position[1], position[2]])
        
        self.unwrapped_positions = np.array(self.unwrapped_positions)
        
    def calculate_center_of_mass(self) -> None:
        """
        Calculate the center of mass of the cluser
        
        Returns:
            None.
        """        
        self.center_of_mass = np.mean(self.unwrapped_positions, axis=0)
    
    def write_coordinates(self, path_to_directory) -> None:
        """
        Write the cluster coordinates to an XYZ file.
        
        Parameters:
            path_to_directory (str): Path to the directory where the XYZ file will be saved.
        
        Returns:
            None.
        """
        if self.size <= 1:
            # Do not write file if size is 1 or inferior
            return
        
        simulation_box = self.box.get_box_dimensions(self.frame)
        
        if not os.path.exists(os.path.join(path_to_directory, "unwrapped_clusters")):
            # Create the output directory if it doesn't exists
            os.makedirs(os.path.join(path_to_directory, "unwrapped_clusters"))
        
        path_to_directory = os.path.join(path_to_directory, "unwrapped_clusters")
        
        filepath = os.path.join(path_to_directory, f"{self.connectivity}_{self.frame}_{self.root_id}.xyz")
        
        with open(filepath, 'w') as output:
            output.write(f"{len(self.atoms)}\nLattice=\"{simulation_box[0]} 0.0 0.0 0.0 {simulation_box[1]} 0.0 0.0 0.0 {simulation_box[2]}\"\n")
            max_length = max([len(str(atom.id)) for atom in self.atoms])
            for atom, pos in zip(self.atoms, self.unwrapped_positions):
                output.write(f"{atom.element} {str(atom.id).ljust(max_length)} {pos[0]:5.5f} {pos[1]:5.5f} {pos[2]:5.5f}\n")

        output.close()
    
    def calculate_gyration_radius(self) -> None:
        """
        Calculate the gyration radius of the cluster.
        
        Returns:
            None.
        """
        self.gyration_radius = 0
        for i in range(self.unwrapped_positions.shape[0]):
            squared_rij = np.linalg.norm(self.unwrapped_positions[i, :] - self.unwrapped_positions[:, :], axis=1)** 2
            self.gyration_radius += np.sum(squared_rij)
        
        # Normalize the sum by 0.5 s²
        self.gyration_radius = np.sqrt((0.5 / (self.size**2)) * self.gyration_radius) 
        
    def calculate_percolation_probability(self) -> None:
        """
        Calculate the percolation probability of the cluster.
        
        Returns:
            None.
        """
        percolation_x = False
        percolation_y = False
        percolation_z = False
        
        box_size = self.box.get_box_dimensions(self.frame)
        
        for i in range(self.unwrapped_positions.shape[0]):
            dx = np.abs(self.unwrapped_positions[i, 0] - self.unwrapped_positions[:, 0])
            dy = np.abs(self.unwrapped_positions[i, 1] - self.unwrapped_positions[:, 1])
            dz = np.abs(self.unwrapped_positions[i, 2] - self.unwrapped_positions[:, 2])
            
            dx = np.max(dx)
            dy = np.max(dy)
            dz = np.max(dz)
            
            if dx > box_size[0]:
                percolation_x = True
            if dy > box_size[1]:
                percolation_y = True
            if dz > box_size[2]:
                percolation_z = True
            
        if percolation_x:
            self.percolation_probability += 'x'
        if percolation_y:
            self.percolation_probability += 'y'
        if percolation_z:
            self.percolation_probability += 'z'
    
    def calculate_order_parameter(self) -> None:
        """
        Calculate the order parameter of the cluster.
        
        Returns:
            None.
        """
        
        if len(self.percolation_probability) == 0:
            return
        elif len(self.percolation_probability) == 1:
            self.order_parameter[0] = self.size / self.number_of_nodes
            self.order_parameter[1] = 0.0
            self.order_parameter[2] = 0.0
        elif len(self.percolation_probability) == 2:
            self.order_parameter[0] = self.size / self.number_of_nodes
            self.order_parameter[1] = self.size / self.number_of_nodes
            self.order_parameter[2] = 0.0
        elif len(self.percolation_probability) == 3:
            self.order_parameter[0] = self.size / self.number_of_nodes
            self.order_parameter[1] = self.size / self.number_of_nodes
            self.order_parameter[2] = self.size / self.number_of_nodes
            
    def calculate_unwrapped_positions(self, criteria, chain, quiet=False) -> None:
        """
        Calculate the unwrapped positions of the atoms in the cluster to make it continuous.
        
        Parameters:
            criteria (str): Criteria to find the clusters ("bond" or "distance").
            chain (list): List of the successive elements forming the clusters.
        
        Returns:
            None.
        """
        # Initialize the dictionary of the unwrapped positions
        stack = [self.atoms[0].parent]
        
        dict_positions = {stack[0].id: stack[0].position}
    
        # Perform DFS traversal to compute unwrapped positions
        box_size = self.box.get_box_dimensions(self.frame)
        
        if criteria == "bond" and len(chain) == 3:
            pass
        elif criteria == "distance" and len(chain) == 2:
            pass
        else:
            raise ValueError(f"\t\tERROR: Something wrong with criteria - chain {criteria}-{chain}.")
        
        if criteria == "bond":
            node_1 = chain[0]
            bridge = chain[1]
            node_2 = chain[2]
            if not quiet:
                while tqdm(stack, desc="Unwrapping clusters ...", colour="YELLOW", leave=False, unit='atom'):
                    current_atom = stack.pop()
                    if current_atom.element == node_1:
                        for first_neighbour in current_atom.neighbours:
                            if first_neighbour.element == bridge:
                                for second_neighbour in first_neighbour.neighbours:
                                    if (
                                        second_neighbour.element == node_2 
                                        and second_neighbour.id not in dict_positions
                                        and second_neighbour.cluster_id == self.root_id
                                    ):
                                        # Compute relative position from the current atom to its second_neighbour
                                        relative_position = self.unwrap_position(
                                            second_neighbour.position - current_atom.position, box_size
                                        )

                                        # Accumulate relative position to get unwrapped position
                                        dict_positions[second_neighbour.id] = (
                                            dict_positions[current_atom.id] + relative_position
                                        )

                                        stack.append(second_neighbour)
            else:
                while stack:
                    current_atom = stack.pop()
                    if current_atom.element == node_1:
                        for first_neighbour in current_atom.neighbours:
                            if first_neighbour.element == bridge:
                                for second_neighbour in first_neighbour.neighbours:
                                    if (
                                        second_neighbour.element == node_2 
                                        and second_neighbour.id not in dict_positions
                                        and second_neighbour.cluster_id == self.root_id
                                    ):
                                        # Compute relative position from the current atom to its second_neighbour
                                        relative_position = self.unwrap_position(
                                            second_neighbour.position - current_atom.position, box_size
                                        )

                                        # Accumulate relative position to get unwrapped position
                                        dict_positions[second_neighbour.id] = (
                                            dict_positions[current_atom.id] + relative_position
                                        )

                                        stack.append(second_neighbour)                            
                            
        if criteria == "distance":
            node_1 = chain[0]
            node_2 = chain[1]
            
            if not quiet:
                while tqdm(stack, desc="Unwrappring clusters ...", colour="YELLOW", leave=False, unit='atom'):
                    current_atom = stack.pop()
                    if current_atom.element == node_1:
                        for neighbour in current_atom.neighbours:
                            if (neighbour.element == node_2
                                and neighbour.id not in dict_positions
                                and neighbour.cluster_id == self.root_id):
                                
                                # Compute relative position from the current atom to its neighbour
                                relative_position = self.unwrap_position(
                                    neighbour.position - current_atom.position, box_size
                                )
                                
                                # Accumulate relative position to get unwrapped position
                                dict_positions[neighbour.id] = (
                                    dict_positions[current_atom.id] + relative_position
                                )
                                
                                stack.append(neighbour)
            else:
                while stack:
                    current_atom = stack.pop()
                    if current_atom.element == node_1:
                        for neighbour in current_atom.neighbours:
                            if (neighbour.element == node_2
                                and neighbour.id not in dict_positions
                                and neighbour.cluster_id == self.root_id):
                                
                                # Compute relative position from the current atom to its neighbour
                                relative_position = self.unwrap_position(
                                    neighbour.position - current_atom.position, box_size
                                )
                                
                                # Accumulate relative position to get unwrapped position
                                dict_positions[neighbour.id] = (
                                    dict_positions[current_atom.id] + relative_position
                                )
                                
                                stack.append(neighbour)
                                
        self.set_indices_and_positions(dict_positions)

    def unwrap_position(self, vector, box_size):
        """
        Unwraps position considering periodic boundary conditions.

        Parameters:
            vector (list): Vector defined the difference of the composants of the vectors position of an atom and its nearest neighbour in the cluster.
            box_size (list): Dimensions of the periodic box in 3D space.

        Returns:
            Unwrapped position as a tuple.
        """
        unwrapped_position = []
        for i in range(3):  # Assuming 3D space
            delta = vector[i] - round(vector[i] / box_size[i]) * box_size[i]
            unwrapped_position.append(delta)
        return tuple(unwrapped_position)