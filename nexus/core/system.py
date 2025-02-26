"""
System module for the Nexus project.

This module defines the System class, which represents a system of atoms and provides methods for analyzing and manipulating the system.

Imports:
    - Standard libraries
    - Third-party libraries
    - Internal modules

Classes:
    - System: Represents a system of atoms and provides methods for analyzing and manipulating the system.
"""

# Standard library imports
import importlib
import os
import re
import inspect

# Third-party imports
import numpy as np
from tqdm import tqdm
from scipy.spatial import cKDTree

# Internal imports
from .cutoff import Cutoff
from .cluster import Cluster
from ..utils.generate_color_gradient import generate_color_gradient

class System:
    """
    Represents a system of atoms and provides methods for analyzing and manipulating the system.

    Attributes:
        settings (Settings): Settings object containing the list of all the parameters.
        atoms (list): List of all the atoms in the system.
        box (Box): The Box object containing the lattice information at each frame.
        clusters (list): List of all the clusters of the system.
        counter_c (int): Counter of Cluster object.
        frame (int): Frame of the system in the trajectory.
        cutoffs (Cutoff): Cutoff object managing cutoff distances for pairs of elements.

    Methods:
        __init__: Initializes a System object.
        add_atom: Adds an Atom object to the list of atoms.
        add_cluster: Adds a Cluster object to the list of clusters.
        get_atoms: Returns the list of atoms.
        get_positions: Returns the list of positions and elements of all Atom objects.
        get_positions_by_element: Returns the list of positions of all Atom objects of the same element.
        get_atoms_by_element: Returns the list of Atom objects belonging to the same species.
        get_unique_element: Returns the unique elements present in the system along with their counts.
        reset_cluster_indexes: Resets the cluster indexes for all Atom objects in the system.
        wrap_atomic_positions: Wraps atomic positions inside the simulation box using periodic boundary conditions.
        compute_mass: Returns the mass of the system in atomic unit.
        calculate_neighbours: Calculates the nearest neighbours of all atoms in the system.
        calculate_structural_units: Determines the structural units and other structural properties.
        get_all_clusters: Returns the list of all Cluster objects associated with the given connectivity.
        get_filtered_clusters: Returns the list of Cluster objects associated with the given connectivity.
        get_all_cluster_sizes: Returns the list of cluster sizes associated with the given connectivity.
        get_filtered_cluster_sizes: Returns the list of cluster sizes associated with the given connectivity.
        write_coordinates_all_in_one: Writes the cluster coordinates to an XYZ file.
        decrypt_connectivity: Decrypts the connectivity string to get the atomic species and the number of neighbors.
        find: Finds the root of the cluster to which the given atom belongs.
        union: Unions the two clusters to which the given atoms belong.
        find_clusters: Finds clusters based on specified criteria.
    """

    def __init__(self, settings) -> None:
        """
        Initializes a System object.

        Parameters:
            settings (Settings): Settings object containing the list of all the parameters.
        """
        self.settings : object = settings   # Settings object containing the list of all the parameters
        self.atoms : list = []              # List of all the atoms 
        self.box : object = None               # The Box object containing the lattice information at each frame
        self.clusters : list = []           # List of the all the clusters of the system
        self.counter_c : int = 0            # Counter of Cluster object
        self.frame : int = 0                # Frame of the system in the trajectory
        self.number_of_nodes : dict = {}    # Number of nodes in all clusters of the system
        
        # Set the cutoffs of the system.
        self.cutoffs : object = Cutoff(settings.cutoffs.get_value()) # Cutoffs of the system
        
    def add_atom(self, atom) -> None:
        """
        Add an Atom object to the list of atoms.
        
        Returns:
            None.
        """
        module = importlib.import_module(f"nexus.extensions.{self.settings.extension.get_value()}")
        transformed_atom = module.transform_into_subclass(atom)
        self.atoms.append(transformed_atom)
    
    def add_cluster(self, cluster:object) -> None:
        """
        Add a Cluster object to the list of clusters.
        
        Returns:
            None.
        """
        self.clusters.append(cluster)
        
    def get_atoms(self) -> list:
        """
        Return the list of atoms.
        
        Returns:
            list : list of Atom objects in the system.
        """
        return self.atoms
    
    def get_positions(self) -> tuple:
        """
        Return the list of positions and elements of all Atom objects.
        
        Returns:
            tuple : the filtered position in a np.array and their associated elements in a np.array.
        """
        filtered_positions = list(
                map(
                    lambda atom: atom.position,
                    filter(
                        lambda atom: hasattr(atom, "frame")
                        and atom.frame == self.frame,
                        self.atoms,
                    ),
                )
            )
        
        filtered_elements = list(
                map(
                    lambda atom: atom.element,
                    filter(
                        lambda atom: hasattr(atom, "frame")
                        and atom.frame == self.frame,
                        self.atoms,
                    ),
                )
            )
        
        return np.array(filtered_positions), np.array(filtered_elements)
        
    def get_positions_by_element(self, element) -> np.array:
        """
        Return the list of positions of all Atom objects of the same element.
        
        Returns:
            np.array : Filtered positions.
        """
        filtered_positions = list(
                map(
                    lambda atom: atom.position,
                    filter(
                        lambda atom: hasattr(atom, "frame")
                        and atom.frame == self.frame
                        and atom.element == element,
                        self.atoms,
                    ),
                )
            )
        
        return np.array(filtered_positions)
    
    def get_atoms_by_element(self, element) -> list:
        """
        Return the list of Atom objects belonging to the same species.
        
        Returns:
            list : list of Atom objects.
        """
        filtered_atoms = list(
                filter(
                    lambda atom: hasattr(atom, "frame", "element")
                    and atom.frame == self.frame
                    and atom.element == element,
                    self.atoms,
                )
            )
        
        return filtered_atoms
    
    def get_unique_element(self) -> np.array:
        """
        Return the uniques elements present in the system along with their counts.
        
        Returns:
            np.array : array of the unique element in the system.
        """
        filtered_elements = np.array(
                list(
                    map(
                        lambda atom: atom.element,
                        filter(
                            lambda atom: hasattr(atom, "frame")
                            and atom.frame == self.frame,
                            self.atoms,
                        ),
                    )
                )
            )
        return np.unique(filtered_elements, return_counts=True)
    
    def reset_cluster_indexes(self) -> None:
        """
        Reset the cluster indexes for all Atom objects in the system.
        
        Returns:
            None.
        """
        for atom in self.atoms:
            atom.reset_cluster()
            
    def wrap_atomic_positions(self) -> None:
        """
        Wrap atomic positions inside the simulation box using the periodic boundary conditions.
        
        Returns:
            None.
        """
        color_gradient = generate_color_gradient(len(self.atoms))
        if self.settings.quiet.get_value() == False:
            progress_bar = tqdm(self.atoms, desc="Wrapping positions inside the box ...", colour="#0dff00", leave=False, unit="atom")
        else:
            progress_bar = self.atoms
        color = 0
        for atom in progress_bar:
            # Updating progress bar
            if self.settings.quiet.get_value() == False:
                progress_bar.set_description(f"Wrapping positions inside the box {atom.id} ...")
                progress_bar.colour = "#%02x%02x%02x" % color_gradient[color]
                color += 1
            
            # Getting box dimensions at the current frame
            box_size = self.box.get_box_dimensions(self.frame)
            
            # Loop over the dimension of the simulation box (ie 3D)
            for i in range(len(box_size)):
                # Apply periodic boundary conditions for each dimension
                atom.position[i] = np.mod(atom.position[i] + box_size[i], box_size[i])
                
    def compute_mass(self) -> float:
        """
        Return the mass of the system in atomic unit.
        
        Returns:
            float : Total mass of the system.
        """
        mass = 0
        for atom in self.atoms:
            mass += atom.atomic_mass
            
        return mass
    
    def calculate_neighbours(self) -> None:
        """
        Calculate the nearest neighbours of all the atom in the system.        
        - NOTE: this method is extension dependant.
        
        Returns:
            None.
        """
        
        # Wrap all the positions inside the simulation box first
        self.wrap_atomic_positions()
        
        # Get the simulation box size
        box_size = self.box.get_box_dimensions(self.frame)
        
        # Get all the atomic positions
        positions, mask = self.get_positions()
        
        # Get the maximum value of the cutoffs of the system
        max_cutoff = self.cutoffs.get_max_cutoff()
        
        # Calculate the tree with the pbc applied
        tree_with_pbc = cKDTree(positions, boxsize=box_size)
        
        # Set the progress bar
        color_gradient = generate_color_gradient(len(positions))
        if self.settings.quiet.get_value() == False:
            progress_bar = tqdm(range(len(positions)), desc="Fetching nearest neighbours ...", colour="#00ffff", leave=False, unit="atom")
        else:
            progress_bar = range(len(positions))
        
        # Loop over the atomic positions
        for i in progress_bar:
            # Update progress bar
            if self.settings.quiet.get_value() == False:
                progress_bar.set_description(f"Fetching nearest neighbours {i} ...")
                progress_bar.colour = "#%02x%02x%02x" % color_gradient[i]
            
            # Process with pbc applied
            # Query the neighbouring atoms within the cutoff distance
            index = tree_with_pbc.query_ball_point(positions[i], max_cutoff)
            
            # Calculate the distance with k nearest neighbours
            distances, indices = tree_with_pbc.query(positions[i], k=len(index))
            
            # Check if result is a list or a int
            if isinstance(indices, int):
                # indices is an int, turn indices into a list of a single int
                indices = [indices]
            
            # Check if results is a list of a int
            if isinstance(distances, int):
                # distances is an int, turn distances into a list of a single int
                distances = [distances]
            
            # Add the nearest neighbours to central atom
            for j in indices:
                self.atoms[i].add_neighbour(self.atoms[j])
            
            self.atoms[i].filter_neighbours(distances)
            self.atoms[i].calculate_coordination()
    
    def calculate_concentrations(self, extension) -> None:
        """
        Determine the structural units and other structural properties.
        - NOTE: this method is extension dependant.
        
        Parameters:
            extension (str) : name of the extension to use to calculate the structural units.
        
        Returns:
            None.
        """
        
        module = importlib.import_module(f"nexus.extensions.{extension}")
        criteria = self.settings.cluster_settings.get_value()['criteria']
        self.concentrations = module.calculate_concentrations(self.get_atoms(), criteria, self.settings.quiet.get_value())
    
    
    #____________CLUSTERS METHODS____________
    
    def set_concentrations(self, connectivity: str) -> None:
        """
        Set the concentrations of the structural units.
        
        Parameters:
            dict_units (dict) : Dict of the structural units.
            cluster_settings (dict) : Dict of the cluster settings.
        
        Returns:
            None.
        """
        
        for connec, concentration in self.concentrations.items():
            list_clusters = self.get_all_clusters(connec)
            for cluster in list_clusters:
                cluster.concentration = concentration

    def get_concentration(self, connectivity: str) -> float:
        """
        Return the concentration of the sites for a given connectivity.
        
        Parameters:
            connectivity (str) : Connectivity of the cluster.
        
        Returns:
            float : Concentration of the sites.
        """
        
        return self.concentrations[connectivity]
                        
    def get_all_clusters(self, connectivity: str) -> list:
        """
        Return the list of all Cluster objects associated with the given connectivity.
        
        Parameters:
            connectivity (str) : Connectivity of the cluster.
            
        Returns:
            list : List of the clusters associated with the given connectivity.
        """
        return [cluster for cluster in self.clusters if cluster.connectivity == connectivity]
    
    def get_filtered_clusters(self, connectivity:str) -> list:
        """
        Return the list of Cluster objects associated with the given connectivity.
        This method excluds clusters of size 1 and percolating (1D, 2D or 3D) clusters.
        
        Parameters:
            connectivity (str) : Connectivity of the cluster.
            
        Returns:
            list : List of the clusters associated with the given connectivity.
        """
        return [cluster for cluster in self.clusters if (
                    (cluster.connectivity == connectivity) 
                and (cluster.size > 1) 
                and (len(cluster.percolation_probability) == 0)
                )
            ]

    def get_all_cluster_sizes(self, connectivity:str) -> list:
        """
        Return the list of cluster size associated with the given connectivity.
        This method excluds clusters of size 1.
        
        Parameters:
            connectivity (str) : Connectivity of the cluster.
            
        Returns:
            list : List of the clusters associated with the given connectivity.
        """
        sizes = [cluster.size for cluster in self.clusters if (
                    (cluster.connectivity == connectivity) 
                and (cluster.size > 1) 
                )
            ]
        if len(sizes) == 0:
            # list of sizes is empty, return [0]
            return [0]
        
        return sizes
    
    def get_filtered_cluster_sizes(self, connectivity:str) -> list:
        """
        Return the list of cluster size associated with the given connectivity.
        This method excluds clusters of size 1 and percolating (1D, 2D or 3D) clusters.
        
        Parameters:
            connectivity (str) : Connectivity of the cluster.
            
        Returns:
            list : List of the clusters associated with the given connectivity.
        """
        sizes = [cluster.size for cluster in self.clusters if (
                    (cluster.connectivity == connectivity) 
                and (cluster.size > 1) 
                and (len(cluster.percolation_probability) == 0)
                )
            ]
        if len(sizes) == 0:
            # list of sizes is empty, return [0]
            return [0]
        
        return sizes
    
    def get_cluster_sizes_distribution(self, connectivity:str) -> dict:
        """
        Return the cluster size distribution of the clusters associated with the given connectivity.
        
        Parameters:
            connectivity (str) : Connectivity of the cluster.
            
        Returns:
            dict : Dict of the cluster size distribution, where keys are the size of the clusters associated with the given connectivity.
        """
        sizes = self.get_filtered_cluster_sizes(connectivity)
        
        dict_sizes = {}
        
        for size in sizes:
            dict_sizes[size] = len([cluster for cluster in self.clusters if (cluster.size == size and cluster.connectivity == connectivity)])
        
        return dict_sizes

    def get_gyration_radius_distribution(self, connectivity:str, list_sizes:list) -> dict:
        """
        Return the gyration radius distribution of the clusters associated with the given connectivity.
        
        Parameters:
            list_sizes (list) : List of the cluster sizes.
            
        Returns:
            dict : Dict of the gyration radius distribution, where keys are the size of the clusters associated with the given connectivity.
        """
        dict_rgyr = {}
        
        for size in list_sizes:
            dict_rgyr[size] = [cluster.gyration_radius for cluster in self.clusters if (cluster.size == size and cluster.connectivity == connectivity)]
        
        return dict_rgyr
    
    def calculate_order_parameter(self, connectivity:str) -> list:
        """
        Calculate the order parameter of the percolating cluster.
        
        Parameters:
            connectivity (str) : Connectivity of the cluster.
        
        Returns:
            list : Order parameter of the percolating cluster [1d value, 2d value, 3d value].
        """
        
        # get the percolating cluster
        cluster = [cluster for cluster in self.clusters if (cluster.connectivity == connectivity) and (len(cluster.percolation_probability) > 0)]
        
        if len(cluster) == 0:
            # No percolating cluster
            return [0.0] * 3
        if len(cluster) > 1:
            # More than one percolating cluster
            # print(f"\tWARNING: More than one percolating cluster found for the connectivity {connectivity}. Taking the biggest one.")
            
            # get the biggest cluster
            cluster = max(cluster, key=lambda x: x.size)
            
            order_parameter = cluster.order_parameter
            
        else:
            order_parameter = cluster[0].order_parameter
        
        return order_parameter
    
    def calculate_percolation_probability(self, connectivity:str) -> list:
        """
        Calculate the percolation probability of the percolating cluster.
        
        Parameters:
            connectivity (str) : Connectivity of the cluster.
            
        Returns:
            list : Percolation probability of the percolating cluster [1d value, 2d value, 3d value].
        """
        
        cluster = [cluster for cluster in self.clusters if (cluster.connectivity == connectivity) and (len(cluster.percolation_probability) > 0)]
        
        if len(cluster) == 0:
            # No percolating cluster
            return [0.0] * 3
        if len(cluster) > 1:
            # More than one percolating cluster
            # print(f"\tWARNING: More than one percolating cluster found for the connectivity {connectivity}. Taking the biggest one.")
            
            # get the biggest cluster
            cluster = max(cluster, key=lambda x: x.size)
            
            percolation_probability = cluster.percolation_probability
            if len(percolation_probability) == 1:
                percolation_probability = [1.0, 0.0, 0.0]
                return percolation_probability
            if len(percolation_probability) == 2:
                percolation_probability = [1.0, 1.0, 0.0]
                return percolation_probability
            if len(percolation_probability) == 3:
                percolation_probability = [1.0, 1.0, 1.0]
                return percolation_probability
        else:
            percolation_probability = cluster[0].percolation_probability
            if len(percolation_probability) == 1:
                percolation_probability = [1.0, 0.0, 0.0]
                return percolation_probability
            if len(percolation_probability) == 2:
                percolation_probability = [1.0, 1.0, 0.0]
                return percolation_probability
            if len(percolation_probability) == 3:
                percolation_probability = [1.0, 1.0, 1.0]
                return percolation_probability

    def write_coordinates_all_in_one(self, connectivity, path_to_directory) -> None:
        """
        Write the cluster coordinates to an XYZ file and a list.txt file with their paths.
        
        Parameters:
            path_to_directory (str): Path to the directory where the XYZ file will be saved.
            
        Returns:
            None.
        """
        
        clusters = [cluster for cluster in self.clusters if cluster.connectivity ==  connectivity and cluster.size > 1]
        
        if len(clusters) == 0:
            #No cluster 
            return
        
        sizes = sum([cluster.size for cluster in clusters if cluster.connectivity ==  connectivity])
        
        indices = ([(max(cluster.indices)) for cluster in clusters if cluster.connectivity == connectivity])
        
        
        simulation_box = self.box.get_box_dimensions(self.frame)
        
        if not os.path.exists(os.path.join(path_to_directory, "unwrapped_clusters")):
            # Create the output directory if it doesn't exists
            os.makedirs(os.path.join(path_to_directory, "unwrapped_clusters"))
        
        path_to_directory = os.path.join(path_to_directory, "unwrapped_clusters")
        
        max_length = len(str(max(indices)))
        
        with open(os.path.join(path_to_directory, f"{connectivity}_{self.frame}_all-in-one.xyz"), 'w') as output:
            output.write(f"{sizes}\nLattice=\"{simulation_box[0]} 0.0 0.0 0.0 {simulation_box[1]} 0.0 0.0 0.0 {simulation_box[2]}\"\n")
        output.close()
        
        with open(os.path.join(path_to_directory, f"{connectivity}_{self.frame}_all-in-one.xyz"), 'a') as output:
            for cluster in clusters:
                for atom, pos in zip(cluster.atoms, cluster.unwrapped_positions):
                    output.write(f"{atom.element} {str(atom.id).ljust(max_length)} {pos[0]:5.5f} {pos[1]:5.5f} {pos[2]:5.5f}\n")
        output.close()
        
    def decrypt_connectivity(self, connectivity):
        """
        Decrypt the connectivity string to get the atomic species and the number of neighbors.
        
        Parameters:
            connectivity (str): The connectivity string to decrypt.
        
        Returns:
            tuple: Tuple containing the atomic species and the number of neighbors.
        """
        
        species = []
        coordinances = re.findall(r'\d+', connectivity)
        units = connectivity.split('-')
        for unit in units:
            matches = re.findall(r'[A-Z][a-z]?', unit)
            for match in matches:
                if len(match) == 2 and match[1].isupper():
                    species.extend(list(match))
                else:
                    species.append(match)
        return np.unique(np.array(species)), int(coordinances[0]), int(coordinances[1])
    
    def find(self, atom) -> object:
        """
        Find the root of the cluster to which the given atom belongs.
        
        Parameters:
            atom (Atom): Atom object for which to find the root.
        
        Returns:
            root (Atom): Root Atom object of the cluster
        """
        if atom.parent != atom:
            atom.parent = self.find(atom.parent)
        return atom.parent

    def union(self, atom_1, atom_2) -> None:
        """
        Union the two clusters to which the given atoms belong.
        
        Parameters:
            atom_1 (Atom): First Atom object.
            atom_2 (Atom): Second Atom object
            
        Returns:
            None.
        """
        root_1 = self.find(atom_1)
        root_2 = self.find(atom_2)
        
        if root_1 != root_2:
            root_2.parent = root_1

    def find_clusters(self, connectivity) -> None:
        """
        Find clusters based on specified criteria.
        
        Parameters:
        - connectivity (str): The connectivity type for which to find the clusters. (e.g "SiO5-SiO5" in the case of SiO2)
        
        Returns:
            None.
        """
        
        # Load cluster properties to analyse
        cluster_settings = self.settings.cluster_settings.get_value()
        criteria = cluster_settings['criteria']
        
        list_of_elements, z1, z2 = self.decrypt_connectivity(connectivity)
        
        if criteria != "distance" and criteria != "bond":
            raise ValueError(f"\tERROR: Criteria '{criteria}' is not supported. Please select 'bond' or 'distance'.")
        for e in list_of_elements:
            if e not in cluster_settings['connectivity']:
                raise ValueError(f"\tERROR: Something wrong in the connectivity in cluster settings.")


        if criteria == "bond":
            node_1  = cluster_settings['connectivity'][0]
            bridge  = cluster_settings['connectivity'][1]
            node_2  = cluster_settings['connectivity'][2]
            chain = [node_1, bridge, node_2] 

        if criteria == "distance":
            node_1 = cluster_settings['connectivity'][0]
            node_2 = cluster_settings['connectivity'][1]
            chain = [node_1, node_2]
        
        networking_atoms = [atom for atom in self.atoms if (atom.element == node_1 or atom.element == node_2) and (atom.coordination == z1 or atom.coordination == z2)]

        number_of_nodes = 0
        
        color_gradient = generate_color_gradient(len(networking_atoms))
        if self.settings.quiet.get_value() == False:
            progress_bar = tqdm(networking_atoms, desc=f"Finding clusters {connectivity} ...", colour="BLUE", leave=False)
        else:
            progress_bar = networking_atoms
        colour = 0
        for atom in progress_bar:
            # Update progress bar
            if self.settings.quiet.get_value() == False:
                progress_bar.set_description(f"Finding clusters {connectivity} ...")
                progress_bar.colour = "#%02x%02x%02x" % color_gradient[colour]
                colour += 1
            
            if criteria == 'distance':
                for neighbour in atom.neighbours:
                    if (atom.element == node_1 and neighbour.element == node_2) and (atom.coordination == z1 and neighbour.coordination == z2):
                        self.union(neighbour, atom)
            if criteria == 'bond':
                for neighbour in atom.neighbours: 
                    if neighbour.element == bridge:
                        for second_neighbour in neighbour.neighbours:
                            if (atom.element == node_1 and second_neighbour.element == node_2) and (atom.coordination == z1 and second_neighbour.coordination == z2):
                                self.union(second_neighbour, atom)
        
        clusters_found = {}
        local_clusters = []
        
        for atom in networking_atoms:
            root = self.find(atom)
            clusters_found.setdefault(root.id, []).append(atom)
        
        color_gradient = generate_color_gradient(len(clusters_found))
        if self.settings.quiet.get_value() == False:
            progress_bar = tqdm(range(len(clusters_found)), desc=f"Calculating clusters properties {connectivity} ...", colour="GREEN", leave=False)
        else:
            progress_bar = range(len(clusters_found))
        colour = 0
        for i in progress_bar:
            cluster = list(clusters_found.values())[i]
            
            # Update progress bar
            if self.settings.quiet.get_value() == False:
                progress_bar.set_description(f"Calculating clusters properties {connectivity} ...")
                progress_bar.colour = "#%02x%02x%02x" % color_gradient[colour]
                colour += 1
            
            for atom in cluster:
                root = self.find(atom)
                break
            
            current_cluster = Cluster(box=self.box,connectivity=connectivity, root_id=self.counter_c, frame=root.frame, size=len(cluster))
            
            for atom in cluster:
                atom.set_cluster(self.counter_c, connectivity)
                current_cluster.add_atom(atom)
                if len(cluster) > 1:
                    number_of_nodes += 1
            
            current_cluster.calculate_unwrapped_positions(criteria, chain, self.settings.quiet.get_value())
            current_cluster.calculate_center_of_mass()
            current_cluster.calculate_gyration_radius()
            current_cluster.calculate_percolation_probability()
            
            if self.settings.print_clusters_positions.get_value():
                current_cluster.write_coordinates(self.settings._output_directory)
            
            self.clusters.append(current_cluster)
            local_clusters.append(current_cluster)
            self.counter_c += 1
            
            for atom in cluster:
                atom.reset_cluster()
            
        if number_of_nodes == 0:
            number_of_nodes = 1 # avoid zero division
        
        for cluster in local_clusters:
            cluster.number_of_nodes = number_of_nodes
            cluster.calculate_order_parameter()
            
        if self.settings.print_clusters_positions.get_value():
            self.write_coordinates_all_in_one(connectivity, self.settings._output_directory)
            
            
    def find_extra_clusters(self) -> None:
        """
        """
        module = importlib.import_module(f"nexus.extensions.{self.settings.extension.get_value()}")
        
        function_list = [o for o in inspect.getmembers(module) if inspect.isfunction(o[1])]
        function_names = [name for name, _ in function_list]
        
        if 'find_extra_clusters' in function_names:
            cluster_settings = self.settings.cluster_settings.get_value()
            criteria = cluster_settings['criteria']
            
            list_extra_clusters, counter_c = module.find_extra_clusters(atoms=self.atoms, box=self.box, counter_c=self.counter_c, settings=self.settings)
            
            
            for cluster in list_extra_clusters:
                self.clusters.append(cluster)
                self.counter_c += 1
                
            if self.settings.print_clusters_positions.get_value():
                for c in module.LIST_OF_EXTRA_CONNECTIVITIES:
                    self.write_coordinates_all_in_one(c, self.settings._output_directory)


