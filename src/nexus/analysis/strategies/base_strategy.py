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

    A clustering strategy defines a specific algorithm for grouping nodes into
    clusters based on a set of criterion. This class provides the fundamental
    interface for all strategies, including a common union-find implementation
    for efficient cluster merging.

    Each concrete strategy must implement the `build_clusters` method, which
    contains the core logic for identifying and forming the clusters.

    Attributes:
        frame (Frame): The simulation frame containing the nodes to be clustered.
        _nodes (List[Node]): A direct reference to the list of nodes in the frame.
        _settings (Settings): The application settings object, providing access to clustering parameters.
    """

    def __init__(self, frame: Frame, settings: Settings) -> None:
        self.frame: Frame = frame
        self._lattice: np.ndarray = self.frame.lattice
        self._nodes: List[Node] = self.frame.nodes
        self._settings: Settings = settings

    def find(self, node: Node) -> Node:
        """
        Finds the root representative of the set containing the given node.
        Implements path compression for optimization.
        """
        if node.parent != node:
            node.parent = self.find(node.parent)
        return node.parent

    def union(self, node_1: Node, node_2: Node) -> None:
        """
        Merges the sets containing two nodes by making one root the parent of the other.
        """
        root_1 = self.find(node_1)
        root_2 = self.find(node_2)

        if root_1 != root_2:
            root_2.parent = root_1

    @abstractmethod
    def build_clusters(self) -> List[Cluster]:
        """
        Executes the clustering algorithm to find and construct all clusters.
        This method must be implemented by all the strategies.
        """
        pass

    def __str__(self) -> str:
        return f"{self.__class__.__name__}({self._settings})"

    def __repr__(self) -> str:
        return self.__str__()
