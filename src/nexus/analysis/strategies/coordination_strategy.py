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
from ...utils.geometry import cartesian_to_fractional
from .base_strategy import BaseClusteringStrategy
from .search.neighbor_searcher import NeighborSearcher


class CoordinationStrategy(BaseClusteringStrategy):
    """
    A clustering strategy that groups nodes based on their coordination number.

    This strategy extends the basic bonding pattern by adding constraints on the 
    coordination number of the networking nodes. A cluster is formed only if the 
    connected nodes meet the specified coordination number in the settings.

    It supports various modes, including :
    - 'pairwise' : pairs with identical coordination number.
    - 'mixing' : pairs with different coordination number.
    - 'alternating' : pairs with alternating coordination number.
    - 'default' : pairs with any coordination number in specified range.
    """
    def __init__(self, frame: Frame, settings: Settings) -> None:
        self.frame: Frame = frame
        self.clusters: List[Cluster] = []
        self._lattice: np.ndarray = self.frame.lattice
        self._nodes: List[Node] = self.frame.nodes
        self._settings: Settings = settings
        self._counter: int = 0
        self._search_mode = "default"
        self._neighbor_searcher = NeighborSearcher(self.frame, self._settings)
        
    def find_neighbors(self) -> None:
        self._neighbor_searcher.execute()

        # Calculate the coordination number
        for i in range(len(self._nodes)):
            self.calculate_coordination(idx=i)

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
        if self._settings.clustering.criteria == 'bond':
            type1 = self._settings.clustering.connectivity[0]
            type2 = self._settings.clustering.connectivity[1]
            type3 = self._settings.clustering.connectivity[2]

            coordination_range = self._settings.clustering.coordination_range
            
            if self._settings.clustering.with_pairwise:
                self._search_mode = "pairwise"
                connectivities = [f"{type1}{type2}_{i}-{type3}{type2}_{i}" for i in range(coordination_range[0], coordination_range[1] + 1)]
            elif self._settings.clustering.with_mixing:
                self._search_mode = "mixing"
                # Generate all possible combinations within the coordination range
                connectivities = []
                for i in range(coordination_range[0], coordination_range[1] + 1):
                    for j in range(i, coordination_range[1] + 1):  # j >= i to avoid duplicates
                        connectivities.append(f"{type1}{type2}_{i}-{type3}{type2}_{j}")
                        if i != j:  # Add the reverse combination if different
                            connectivities.append(f"{type3}{type2}_{j}-{type1}{type2}_{i}")
            elif self._settings.clustering.with_alternating:
                self._search_mode = "alternating"
                connectivities = []
                for i in range(coordination_range[0], coordination_range[1] + 1):
                    connectivities.append(f"{type1}{type2}_{i}-{type3}{type2}_{i}")
                    if i+1 <= coordination_range[1]:
                        connectivities.append(f"{type3}{type2}_{i}-{type1}{type2}_{i+1}")
            else:
                self._search_mode = "default"
                connectivities = [self._settings.clustering.with_connectivity_name]
        else:
            type1 = self._settings.clustering.connectivity[0]
            type2 = self._settings.clustering.connectivity[1]

            coordination_range = self._settings.clustering.coordination_range
            
            if self._settings.clustering.with_pairwise:
                self._search_mode = "pairwise"
                connectivities = [f"{type1}_{i}-{type2}_{i}" for i in range(coordination_range[0], coordination_range[1] + 1)]
            elif self._settings.clustering.with_mixing:
                self._search_mode = "mixing"
                # Generate all possible combinations within the coordination range
                connectivities = []
                for i in range(coordination_range[0], coordination_range[1] + 1):
                    for j in range(i, coordination_range[1] + 1):  # j >= i to avoid duplicates
                        connectivities.append(f"{type1}_{i}-{type2}_{j}")
                        if i != j:  # Add the reverse combination if different
                            connectivities.append(f"{type2}_{j}-{type1}_{i}")
            elif self._settings.clustering.with_alternating:
                self._search_mode = "alternating"
                connectivities = []
                for i in range(coordination_range[0], coordination_range[1] + 1):
                    connectivities.append(f"{type1}_{i}-{type2}_{i}")
                    if i+1 <= coordination_range[1]:
                        connectivities.append(f"{type2}_{i}-{type1}_{i+1}")
            else:
                self._search_mode = "default"
                connectivities = [self._settings.clustering.with_connectivity_name]

        return connectivities

    def build_clusters(self) -> None:
        # Select the networking nodes based on clustering settings
        # 1 - check node types
        if self._settings.clustering.criteria == 'bond':
            networking_nodes = [node for node in self._nodes if node.symbol in self._settings.clustering.node_types and node.symbol != self._settings.clustering.connectivity[1]]
        else:
            networking_nodes = [node for node in self._nodes if node.symbol in self._settings.clustering.node_types]
        
        # 2 - generate connectivities based on coordination number range
        connectivities = self.get_connectivities()
        
        # 3 - generate clusters based on connectivities
        if self._search_mode == "default":
            self._find_cluster(networking_nodes, connectivities[0], 0, 0)
        else:
            for connectivity in connectivities:
                z1, z2 = connectivity.split('-')
                z1 = int(z1.split('_')[1])
                z2 = int(z2.split('_')[1])
                self._find_cluster(networking_nodes, connectivity, z1, z2)
        
        # 4 - return clusters
        return self.clusters

    def _find_cluster(self, networking_nodes: List[Node], connectivity: str, z1: int, z2: int) -> None:
        number_of_nodes = 0

        lbound = self._settings.clustering.coordination_range[0]
        ubound = self._settings.clustering.coordination_range[1]
        if lbound == ubound:
            coordination_range = np.array([lbound])
        else:
            coordination_range = np.arange(lbound, ubound + 1)
        
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
                            if self._search_mode == "default":
                                if (node.symbol == type1 and neighbor2.symbol == type3) and (node.coordination in coordination_range and neighbor2.coordination in coordination_range):
                                    self.union(neighbor2, node)
                            else:
                                if (node.symbol == type1 and neighbor2.symbol == type3) and (node.coordination == z1 and neighbor2.coordination == z2):
                                    self.union(neighbor2, node)
        
        elif self._settings.clustering.criteria == 'distance':
            type1 = self._settings.clustering.connectivity[0]
            type2 = self._settings.clustering.connectivity[1]
            
            for node in progress_bar:
                for neighbor in node.neighbors:
                    if self._search_mode == "default":
                        if (node.symbol == type1 and neighbor.symbol == type2) and (node.coordination in coordination_range and neighbor.coordination in coordination_range):
                            self.union(neighbor, node)
                    else:
                        if (node.symbol == type1 and neighbor.symbol == type2) and (node.coordination == z1 and neighbor.coordination == z2):
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