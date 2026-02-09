from typing import List
import numpy as np
from tqdm import tqdm
import shutil

# Internal imports
from ...core.node import Node
from ...core.cluster import Cluster
from ...core.frame import Frame
from ...config.settings import Settings
from .base_strategy import BaseClusteringStrategy
from .search.neighbor_searcher import NeighborSearcher


class DistanceStrategy(BaseClusteringStrategy):
    """
    Clustering strategy based on a direct distance criterion.

    Connects any two nodes of specified types that are within a cutoff
    distance of each other. Requires a 2-element connectivity specification
    (e.g., ``["Si", "Si"]``).

    Attributes:
        clusters (List[Cluster]): Clusters accumulated across calls.
        _counter (int): Running count of clusters found.
        _neighbor_searcher (NeighborSearcher): KD-tree based neighbor finder.
    """

    def __init__(self, frame: Frame, settings: Settings) -> None:
        """
        Initialize the strategy.

        Args:
            frame (Frame): The simulation frame to operate on.
            settings (Settings): Configuration settings.
        """
        self.frame: Frame = frame
        self.clusters: List[Cluster] = []
        self._lattice: np.ndarray = self.frame.lattice
        self._nodes: List[Node] = self.frame.nodes
        self._settings: Settings = settings
        self._counter: int = 0
        self._neighbor_searcher = NeighborSearcher(self.frame, self._settings)

    def find_neighbors(self) -> None:
        """Populate each node's neighbor list using the KD-tree searcher."""
        self._neighbor_searcher.execute()

    def get_connectivities(self) -> List[str]:
        """
        Build the connectivity label from the 2-element connectivity setting.

        Returns:
            List[str]: Single-element list with the connectivity label.

        Raises:
            ValueError: If connectivity is not a 2-element list.
        """
        connectivity = self._settings.clustering.connectivity
        if isinstance(connectivity, list) and len(connectivity) == 2:
            return [f"{connectivity[0]}-{connectivity[1]}"]
        else:
            raise ValueError(
                "Connectivity for clustering based on distance criterion must be a list of two elements."
            )

    def build_clusters(self) -> List[Cluster]:
        """
        Build clusters of directly connected networking nodes.

        Applies union-find over all networking node pairs within the cutoff
        distance, then computes cluster properties for each cluster.

        Returns:
            List[Cluster]: The clusters found.
        """
        networking_nodes = [
            node
            for node in self._nodes
            if node.symbol in self._settings.clustering.node_types
        ]
        connectivity = self._settings.clustering.connectivity

        number_of_nodes = 0

        progress_bar_kwargs = {
            "disable": not self._settings.verbose,
            "leave": False,
            "ncols": shutil.get_terminal_size().columns,
            "colour": "green",
        }
        progress_bar = tqdm(
            networking_nodes, desc="Finding clusters ...", **progress_bar_kwargs
        )

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

        progress_bar = tqdm(
            range(len(clusters_found)),
            desc="Calculating clusters properties ...",
            **progress_bar_kwargs,
        )

        for i in progress_bar:
            cluster = list(clusters_found.values())[i]

            root = None
            for node in cluster:
                root = self.find(node)
                break
            if root is None:
                raise ValueError("Cluster has no root")

            current_cluster = Cluster(
                connectivity=self.get_connectivities()[0],
                root_id=root.node_id,
                size=len(cluster),
                settings=self._settings,
                lattice=self._lattice,
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
            number_of_nodes = 1  # avoid zero division

        for cluster in local_clusters:
            cluster.total_nodes = number_of_nodes
            cluster.calculate_unwrapped_positions()
            cluster.calculate_center_of_mass()
            cluster.calculate_gyration_radius()
            cluster.calculate_percolation_probability()
            cluster.calculate_concentration()
            cluster.calculate_order_parameter()

        return self.clusters

