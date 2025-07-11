from typing import List
from scipy.spatial import cKDTree
import numpy as np
import os
from tqdm import tqdm

# Internal imports
from ...core.node import Node
from ...core.frame import Frame
from ...core.cluster import Cluster
from ...config.settings import Settings
from .base_finder import BaseFinder
from ...utils.geometry import cartesian_to_fractional


class SharedBasedFinder(BaseFinder):
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

            # Calculate the coordination number
            self.calculate_coordination(idx=i)
            

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

    def calculate_coordination(self, idx: int) -> None:
        node = self._nodes[idx]
        
        mode = self._settings.clustering.coordination_mode

        # "all_types", "same_type", "different_type", "<node_type>"
        if mode == 'all_types':
            node.set_coordination(len(node.neighbors))
        elif mode == 'same_type':
            node.set_coordination(len([n for n in node.neighbors if n.symbol == node.symbol]))
        elif mode == 'different_type':
            node.set_coordination(len([n for n in node.neighbors if n.symbol != node.symbol]))
        else:
            node.set_coordination(len([n for n in node.neighbors if n.symbol == mode]))

    def get_number_of_shared(self, node_1: Node, node_2: Node) -> int:
        mode = self._settings.clustering.shared_mode

        if mode == 'all_types':
            return len([n for n in node_1.neighbors if n in node_2.neighbors])
        elif mode == 'same_type':
            return len([n for n in node_1.neighbors if n.symbol == node_1.symbol and n in node_2.neighbors])
        elif mode == 'different_type':
            return len([n for n in node_1.neighbors if n.symbol != node_1.symbol and n in node_2.neighbors])
        else:
            return len([n for n in node_1.neighbors if n.symbol == mode and n in node_2.neighbors])

    def find(self, node: Node) -> Node:
        if node.parent != node:
            node.parent = self.find(node.parent)
        return node.parent

    def union(self, node_1: Node, node_2: Node) -> None:
        root_1 = self.find(node_1)
        root_2 = self.find(node_2)
        
        if root_1 != root_2:
            root_2.parent = root_1

    def get_connectivities(self) -> List[str]:
        cn = self._settings.clustering.shared_threshold
        if self._settings.clustering.criteria == 'bond':
            type1 = self._settings.clustering.connectivity[0]
            type2 = self._settings.clustering.connectivity[1]
            type3 = self._settings.clustering.connectivity[2]

            coordination_range = self._settings.clustering.coordination_range
            if self._settings.clustering.with_alternating:
                connectivities = []
                for i in range(coordination_range[0], coordination_range[1] + 1):
                    connectivities.append(f"{type1}{type2}_{i}-{type3}{type2}_{i}-{cn}plus_common_neighbors")
                    if i+1 <= coordination_range[1]:
                        connectivities.append(f"{type3}{type2}_{i}-{type1}{type2}_{i+1}-{cn}plus_common_neighbors")
            else:
                connectivities = [f"{type1}{type2}_{i}-{type3}{type2}_{i}-{cn}plus_common_neighbors" for i in range(coordination_range[0], coordination_range[1] + 1)]
        else:
            type1 = self._settings.clustering.connectivity[0]
            type2 = self._settings.clustering.connectivity[1]

            coordination_range = self._settings.clustering.coordination_range            
            if self._settings.clustering.with_alternating:
                connectivities = []
                for i in range(coordination_range[0], coordination_range[1] + 1):
                    connectivities.append(f"{type1}_{i}-{type2}_{i}-{cn}plus_common_neighbors")
                    if i+1 <= coordination_range[1]:
                        connectivities.append(f"{type2}_{i}-{type1}_{i+1}-{cn}plus_common_neighbors")
            else:
                connectivities = []
                for i in range(coordination_range[0], coordination_range[1] + 1):
                    connectivities.append(f"{type1}_{i}-{type2}_{i}-{cn}plus_common_neighbors")

        return connectivities

    def find_clusters(self) -> None:
        # Select the networking nodes based on clustering settings
        # 1 - check node types
        if self._settings.clustering.criteria == 'bond':
            networking_nodes = [node for node in self._nodes if node.symbol in self._settings.clustering.node_types and node.symbol != self._settings.clustering.connectivity[1]]
        else:
            networking_nodes = [node for node in self._nodes if node.symbol in self._settings.clustering.node_types]
        
        # 2 - generate connectivities based on coordination number range
        connectivities = self.get_connectivities()
        
        # 3 - generate clusters based on connectivities
        for connectivity in connectivities:
            z1, z2, dump = connectivity.split('-')
            z1 = int(z1.split('_')[1])
            z2 = int(z2.split('_')[1])
            self._find_cluster(networking_nodes, connectivity, z1, z2)
        
        # 4 - return clusters
        return self.clusters

    def _find_cluster(self, networking_nodes: List[Node], connectivity: str, z1: int, z2: int) -> None:
        number_of_nodes = 0

        progress_bar_kwargs = {
            "disable": not self._settings.verbose,
            "leave": False,
            "ncols": os.get_terminal_size().columns,
            "colour": "blue"
        }
        progress_bar = tqdm(networking_nodes, desc=f"Finding clusters {connectivity} ...", **progress_bar_kwargs)

        if self._settings.clustering.criteria == 'bond':
            type1 = self._settings.clustering.connectivity[0]
            type2 = self._settings.clustering.connectivity[1]
            type3 = self._settings.clustering.connectivity[2]
            
            for node in progress_bar:
                for neighbor in node.neighbors:
                    if neighbor.symbol == type2:
                        for neighbor2 in neighbor.neighbors:
                            if (node.symbol == type1 and neighbor2.symbol == type3) and (node.coordination == z1 and neighbor2.coordination == z2):
                                if self.get_number_of_shared(node, neighbor2) >= self._settings.clustering.shared_threshold:
                                    self.union(neighbor2, node)
        
        elif self._settings.clustering.criteria == 'distance':
            type1 = self._settings.clustering.connectivity[0]
            type2 = self._settings.clustering.connectivity[1]
            
            for node in progress_bar:
                for neighbor in node.neighbors:
                    if (node.symbol == type1 and neighbor.symbol == type2) and (node.coordination == z1 and neighbor.coordination == z2):
                        if self.get_number_of_shared(node, neighbor) >= self._settings.clustering.shared_threshold:
                            self.union(neighbor, node)
        
        clusters_found = {}
        local_clusters = []

        for node in networking_nodes:
            root = self.find(node)
            clusters_found.setdefault(root.node_id, []).append(node)

        progress_bar = tqdm(range(len(clusters_found)), desc=f"Calculating clusters {connectivity} properties ...", **progress_bar_kwargs)  

        for i in progress_bar:
            cluster = list(clusters_found.values())[i]

            for node in cluster:
                root = self.find(node)
                break

            current_cluster = Cluster(
                connectivity=connectivity,
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
            