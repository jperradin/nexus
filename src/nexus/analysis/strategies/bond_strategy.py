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
from ...utils.geometry import cartesian_to_fractional
from .base_strategy import BaseClusteringStrategy
from .search.neighbor_searcher import NeighborSearcher


class BondingStrategy(BaseClusteringStrategy):
    """
    A clustering strategy that connects two nodes via a 'bridge' of a third node.

    This strategy identifies clusters based on a specified bonding pattern, typically
    a three-node connectivity pattern. For example, a connectivity of ['Si', 'O', 'Si']
    will connect silicon atoms that are sharing a common oxygen atom.
    """
    def __init__(self, frame: Frame, settings: Settings) -> None:
        self.frame: Frame = frame
        self.clusters: List[Cluster] = []
        self._lattice: np.ndarray = self.frame.lattice
        self._nodes: List[Node] = self.frame.nodes
        self._settings: Settings = settings
        self._counter: int = 0
        self._neighbor_searcher = NeighborSearcher(self.frame, self._settings)
        
    def find_neighbors(self) -> None:
        self._neighbor_searcher.execute()
        
    def get_connectivities(self) -> List[str]:
        connectivity = self._settings.clustering.connectivity
        if isinstance(connectivity, list) and len(connectivity) == 3:
            connectivity = [f"{connectivity[0]}-{connectivity[1]}-{connectivity[2]}"]
            return connectivity
        else:
            raise ValueError("Connectivity for clustering based on bond criterion must be a list of three elements.")
            
    
    def build_clusters(self) -> List[Cluster]:
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