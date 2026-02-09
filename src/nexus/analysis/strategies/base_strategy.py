from abc import ABC, abstractmethod
from typing import List
import numpy as np

from ...core.frame import Frame
from ...core.node import Node
from ...config.settings import Settings
from ...core.cluster import Cluster


class BaseClusteringStrategy(ABC):
    """
    Abstract base class for all clustering strategies.

    Defines the interface for grouping nodes into clusters using a union-find
    data structure with path compression. Concrete strategies must implement
    ``build_clusters()`` with the specific algorithm for identifying and
    forming clusters.

    Attributes:
        frame (Frame): The simulation frame containing the nodes to cluster.
        _lattice (np.ndarray): Lattice matrix from the frame.
        _nodes (List[Node]): Direct reference to the nodes in the frame.
        _settings (Settings): Configuration settings with clustering parameters.
    """

    def __init__(self, frame: Frame, settings: Settings) -> None:
        """
        Initialize the strategy.

        Args:
            frame (Frame): The simulation frame to operate on.
            settings (Settings): Configuration settings.
        """
        self.frame: Frame = frame
        self._lattice: np.ndarray = self.frame.lattice
        self._nodes: List[Node] = self.frame.nodes
        self._settings: Settings = settings

    def find(self, node: Node) -> Node:
        """
        Find the root representative of the set containing *node*.

        Uses path compression so that subsequent lookups are O(1).

        Args:
            node (Node): The node to look up.

        Returns:
            Node: The root representative of the set.
        """
        if node.parent != node:
            node.parent = self.find(node.parent)
        return node.parent

    def union(self, node_1: Node, node_2: Node) -> None:
        """
        Merge the sets containing two nodes.

        Args:
            node_1 (Node): First node.
            node_2 (Node): Second node.
        """
        root_1 = self.find(node_1)
        root_2 = self.find(node_2)

        if root_1 != root_2:
            root_2.parent = root_1

    @abstractmethod
    def build_clusters(self) -> List[Cluster]:
        """
        Execute the clustering algorithm and return the resulting clusters.

        Returns:
            List[Cluster]: The clusters found by the strategy.
        """
        pass

    def __str__(self) -> str:
        """Return a string representation including settings."""
        return f"{self.__class__.__name__}({self._settings})"

    def __repr__(self) -> str:
        """Return a string representation for debugging."""
        return self.__str__()
