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


class DistanceStrategy(BaseClusteringStrategy):
    """
    A clustering strategy that connects nodes based on a direct distance criterion.

    This strategy forms clusters by connecting any two nodes of specified types
    that are within a given cutoff distance of each other. It is the simplest
    form of clustering, suitable for identifying aggregates, micelles, or
    simple molecular grouping in systems.
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
        if isinstance(connectivity, list) and len(connectivity) == 2:
            return [f"{connectivity[0]}-{connectivity[1]}"]
        else:
            raise ValueError("Connectivity for clustering based on distance criterion must be a list of two elements.")
            
    
    def build_clusters(self) -> List[Cluster]:
        networking_nodes = [node for node in self._nodes if node.symbol in self._settings.clustering.node_types]
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
                        self.union(neighbor, node)
        
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