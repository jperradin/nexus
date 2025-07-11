from typing import List
from scipy.spatial import cKDTree
import numpy as np
import os
from tqdm import tqdm

# Internal imports
from ...core.node import Node
from ...core.cluster import Cluster
from ...core.frame import Frame
from ...config.settings import Settings
from .base_finder import BaseFinder
from ...utils.geometry import cartesian_to_fractional


class GeneralBondFinder(BaseFinder):
    def __init__(self, frame: Frame, settings: Settings) -> None:
        self.frame: Frame = frame
        self.clusters: List[Cluster] = []
        self._lattice: np.ndarray = self.frame.lattice
        self._nodes: List[Node] = self.frame.nodes
        self._settings: Settings = settings
        self._counter: int = 0
        
    def find_neighbors(self) -> None:
        
        # Get the lattice information 
        lattice = self._lattice

        # Get the wrapped positions of all nodes
        positions = self.frame.get_wrapped_positions()

        # Get the cutoffs of the system
        cutoffs = self._settings.clustering.cutoffs

        # Get the maximum value of the cutoffs of the system
        max_cutoff = 0.0
        for cutoff in cutoffs:
            if cutoff.distance > max_cutoff:
                max_cutoff = cutoff.distance

        # Calculate the graph 
        if self._settings.apply_pbc:
            # Convert positions to fractional
            positions = cartesian_to_fractional(positions, lattice)
            # Create the kdtree with pbc inside the unit cell
            kdtree = cKDTree(positions, boxsize=[1, 1, 1])
        else:
            kdtree = cKDTree(positions)

        progress_bar_kwargs = {
            "disable": not self._settings.verbose,
            "leave": False,
            "ncols": os.get_terminal_size().columns,
            "colour": "green"
        }
        progress_bar = tqdm(range(len(positions)), desc="Fetching nearest neighbors ...", **progress_bar_kwargs)

        # Loop over the node positions
        for i in progress_bar:
            # Query the neighboring nodes within the cutoff distance
            if self._settings.apply_pbc:
                # Convert cutoff to fractional
                fractional_cutoff = max_cutoff / np.linalg.norm(lattice, axis=0).max()
                index = kdtree.query_ball_point(positions[i], fractional_cutoff)
            else:
                index = kdtree.query_ball_point(positions[i], max_cutoff)

            # Calculate the distance with k nearest neighbors
            if self._settings.apply_pbc:
                distances, indices = kdtree.query(positions[i], k=len(index))
                # Convert distances to cartesian
                distances *= np.linalg.norm(lattice, axis=0).max()
            else:
                distances, indices = kdtree.query(positions[i], k=len(index))

            # Check if result is a list or a int
            if isinstance(indices, int):
                # indices is an int, turn indices into a list of a single int
                indices = [indices]

            # Check if result is a list or a int
            if isinstance(distances, float):
                # distances is a float, turn distances into a list of a single float
                distances = [distances]

            # Add the distances and indices to the node
            self._nodes[i].distances = distances
            self._nodes[i].indices = indices

            # Add the nearest neighbors to the node
            for j in indices:
                self._nodes[i].add_neighbor(self._nodes[j])
            
            # Filter the neighbors
            self.filter_neighbors(idx=i, distances=distances)

    def filter_neighbors(self, idx: int, distances: List[float]) -> None:
        new_list_neighbors = []
        new_list_distances  = []
        
        node = self._nodes[idx]

        for k, neighbor in enumerate(node.neighbors):
            rcut = self._settings.clustering.get_cutoff(node.symbol, neighbor.symbol)
            
            if isinstance(distances, float):
                # if 'distances' is a float, it means that the neighbor of this node is itself.
                current_distance = distances
            else:
                current_distance = distances[k]
            
            if current_distance > rcut: # neighbor is too far 
                continue # go next neighbor
            elif current_distance == 0: # neighbor is this node.
                continue # go next neighbor
            else:
                if neighbor.node_id == node.node_id:
                    continue # go next neighbor
                new_list_neighbors.append(neighbor) # keep the neighbor
                new_list_distances.append(current_distance)

        node.neighbors = new_list_neighbors
        node.distances = new_list_distances
        node.indices = [node.node_id for node in new_list_neighbors]
        

    def get_connectivities(self) -> List[str]:
        connectivity = self._settings.clustering.connectivity
        if isinstance(connectivity, list) and len(connectivity) == 3:
            connectivity = [f"{connectivity[0]}-{connectivity[1]}-{connectivity[2]}"]
            return connectivity
        else:
            raise ValueError("Connectivity for clustering based on bond criteria must be a list of three elements.")
            
    
    def find_clusters(self) -> List[Cluster]:
        networking_nodes = [node for node in self._nodes if node.symbol in self._settings.clustering.node_types and node.symbol != self._settings.clustering.connectivity[1]]
        connectivity = self._settings.clustering.connectivity

        number_of_nodes = 0
        
        progress_bar_kwargs = {
            "disable": not self._settings.verbose,
            "leave": False,
            "ncols": os.get_terminal_size().columns,
            "colour": "green"
        }
        progress_bar = tqdm(networking_nodes, desc="Finding clusters ...", **progress_bar_kwargs)

        for node in progress_bar:
            if node.symbol == connectivity[0]:
                for neighbor in node.neighbors:
                    if neighbor.symbol == connectivity[1]:
                        for neighbor2 in neighbor.neighbors:
                            if neighbor2.symbol == connectivity[2]:
                                self.union(neighbor2, node)
        
        clusters_found = {}
        local_clusters = []

        for node in networking_nodes:
            root = self.find(node)
            clusters_found.setdefault(root.node_id, []).append(node)

        progress_bar = tqdm(range(len(clusters_found)), desc="Calculating clusters properties ...", **progress_bar_kwargs)  

        for i in progress_bar:
            cluster = list(clusters_found.values())[i]

            for node in cluster:
                root = self.find(node)
                break

            current_cluster = Cluster(
                connectivity=self.get_connectivities()[0],
                root_id=root.node_id,
                size=len(cluster),
                settings=self._settings,
                lattice=self._lattice
            )

            for node in cluster:
                current_cluster.add_node(node)
                if len(cluster) > 1:
                    number_of_nodes += 1
            
            if len(cluster) > 1:
                self.clusters.append(current_cluster)
                local_clusters.append(current_cluster)
                self._counter += 1

            for node in cluster:
                node.reset_parent()
            
        if number_of_nodes == 0:
            number_of_nodes = 1 # avoid zero division
            
        for cluster in local_clusters:
            cluster.total_nodes = number_of_nodes
            cluster.calculate_unwrapped_positions()
            cluster.calculate_center_of_mass()
            cluster.calculate_gyration_radius()
            cluster.calculate_percolation_probability()
            cluster.calculate_concentration()
            cluster.calculate_order_parameter()

        return self.clusters