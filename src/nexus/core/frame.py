import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Optional

from .node import Node
from ..utils.geometry import wrap_positions
from ..core.cluster import Cluster
from ..config.settings import Settings


@dataclass(slots=True)
class Frame:
    """
    Represents a single snapshot of a trajectory.

    A frame holds the raw node data, the simulation cell lattice, and the clusters
    produced by a clustering strategy. It is created by a reader, populated during
    the analysis pipeline, and consumed by analyzers.

    Attributes:
        frame_id (int): Sequential identifier of this frame in the trajectory.
        nodes (List[Node]): Nodes present in this frame after species filtering.
        lattice (np.ndarray): 3x3 lattice matrix defining the simulation cell.
        _data (Dict[str, np.ndarray]): Internal storage for raw node data (symbols and
            positions) as parsed by the reader.
        _settings (Settings): Configuration settings used for node initialization and
            analysis.
        clusters (Optional[List[Cluster]]): Clusters identified in this frame, or None
            if clustering has not been performed yet.
        connectivities (Optional[List[str]]): List of connectivity labels found in this
            frame, or None if not yet determined.
    """

    frame_id: int
    nodes: List[Node]
    lattice: np.ndarray
    _data: Dict[str, np.ndarray]
    _settings: Settings
    clusters: Optional[List[Cluster]] = None
    connectivities: Optional[List[str]] = None

    def __post_init__(self):
        """Validate field types after dataclass initialization."""
        if not isinstance(self.nodes, list):
            raise TypeError("nodes must be a list of Nodes")
        if self.lattice is not None and not isinstance(self.lattice, np.ndarray):
            raise TypeError("lattice must be a numpy array")

    def initialize_nodes(self) -> None:
        """
        Build the node list from raw data, keeping only selected species.

        Reads symbols and positions from the internal ``_data`` dictionary and creates
        a ``Node`` for each entry whose symbol appears in the configured ``node_types``.
        Nodes are assigned sequential IDs starting from zero.
        """
        id = 0
        symbols = self._data["symbol"]
        positions = self._data["position"]

        if len(symbols) != len(positions):
            raise ValueError("symbols and positions must have the same length")

        for symbol, position in zip(symbols, positions):
            if symbol not in self._settings.clustering.node_types:
                continue
            self.nodes.append(Node(node_id=id, symbol=symbol, position=position))
            id += 1

    def set_lattice(self, lattice: np.ndarray) -> None:
        """
        Set the lattice matrix for this frame.

        Args:
            lattice (np.ndarray): A 3x3 non-singular matrix defining the simulation cell.

        Raises:
            ValueError: If the matrix is not 3x3 or is singular.
        """
        if lattice.shape != (3, 3):
            raise ValueError("lattice must be a 3x3 numpy array")

        try:
            np.linalg.inv(lattice)
        except np.linalg.LinAlgError:
            raise ValueError("lattice must be a non-singular matrix")

        self.lattice = lattice

    def get_lattice(self) -> Optional[np.ndarray]:
        """
        Return the lattice matrix of this frame.

        Returns:
            Optional[np.ndarray]: The 3x3 lattice matrix, or None if not set.
        """
        return self.lattice

    def get_unique_elements(self) -> List[str]:
        """
        Return the unique chemical symbols present in this frame.

        Returns:
            List[str]: Sorted list of unique element symbols.
        """
        return np.unique([node.symbol for node in self.nodes])

    def get_node_by_id(self, node_id: int) -> Optional[Node]:
        """
        Look up a node by its identifier.

        Args:
            node_id (int): The identifier of the node to retrieve.

        Returns:
            Optional[Node]: The matching node, or None if no node has the given ID.
        """
        for node in self.nodes:
            if node.node_id == node_id:
                return node
        return None

    def get_positions(self) -> np.ndarray:
        """
        Return the Cartesian positions of all nodes.

        Returns:
            np.ndarray: Position array of shape (N, 3).
        """
        return np.array([node.position for node in self.nodes])

    def get_positions_by_element(self) -> Dict[str, np.ndarray]:
        """
        Return the Cartesian positions grouped by element symbol.

        Returns:
            Dict[str, np.ndarray]: Mapping of element symbols to position arrays.
        """
        return {
            node.symbol: np.array(
                [node.position for node in self.nodes if node.symbol == node.symbol]
            )
            for node in self.nodes
        }

    def get_wrapped_positions(self) -> np.ndarray:
        """
        Return all node positions wrapped into the simulation cell.

        Returns:
            np.ndarray: Wrapped position array of shape (N, 3).
        """
        return wrap_positions(self.get_positions(), self.lattice)

    def get_wrapped_positions_by_element(self) -> Dict[str, np.ndarray]:
        """
        Return wrapped node positions grouped by element symbol.

        Returns:
            Dict[str, np.ndarray]: Mapping of element symbols to wrapped position arrays.
        """
        return {
            node.symbol: wrap_positions(
                np.array(
                    [node.position for node in self.nodes if node.symbol == node.symbol]
                ),
                self.lattice,
            )
            for node in self.nodes
        }

    def get_clusters(self) -> List[Cluster]:
        """
        Return the clusters identified in this frame.

        Returns:
            List[Cluster]: Clusters belonging to this frame.
        """
        return self.clusters

    def get_nodes(self) -> List[Node]:
        """
        Return the list of nodes in this frame.

        Returns:
            List[Node]: Nodes belonging to this frame.
        """
        return self.nodes

    def get_networking_nodes(self) -> int:
        """
        Return the total number of nodes that belong to any cluster.

        Returns:
            int: Sum of all cluster sizes in this frame.
        """
        total_sizes = [c.get_size() for c in self.clusters]
        return np.sum(total_sizes)

    def get_connectivities(self) -> List[str]:
        """
        Return the connectivity labels found in this frame.

        Returns:
            List[str]: Connectivity labels (e.g., ["Si-Si", "Si-O-Si"]).
        """
        return self.connectivities

    def set_connectivities(self, connectivities: List[str]) -> None:
        """
        Set the connectivity labels for this frame.

        Args:
            connectivities (List[str]): Connectivity labels to assign.
        """
        self.connectivities = connectivities

    def add_cluster(self, cluster: Cluster) -> None:
        """
        Append a cluster to this frame, initializing the cluster list if needed.

        Args:
            cluster (Cluster): The cluster to add.
        """
        if self.clusters is None:
            self.clusters = []
        self.clusters.append(cluster)

    def set_clusters(self, clusters: List[Cluster]) -> None:
        """
        Replace the cluster list and update each cluster's lattice.

        Args:
            clusters (List[Cluster]): Clusters to assign to this frame.
        """
        for cluster in clusters:
            cluster.frame_id = self.frame_id
            cluster.set_lattice(self.lattice)
        self.clusters = clusters

    def get_concentration(self) -> Dict[str, float]:
        """
        Compute the concentration of networking nodes for each connectivity.

        For each connectivity label, calculates the ratio of networking nodes to the total
        number of nodes whose symbols match the connectivity endpoints.

        Returns:
            Dict[str, float]: Mapping of connectivity labels to their concentration values.
        """
        concentrations = {}
        for connectivity in self.connectivities:
            for cluster in self.clusters:
                if cluster.get_connectivity() == connectivity:
                    connectivity_settings = self._settings.clustering.connectivity
                    networking_nodes = len(
                        [
                            node
                            for node in self.nodes
                            if node.symbol
                            in [connectivity_settings[0], connectivity_settings[-1]]
                        ]
                    )
                    concentrations[connectivity] = (
                        cluster.total_nodes / networking_nodes
                    )
                    break

        for connectivity in self.connectivities:
            if connectivity not in concentrations:
                concentrations[connectivity] = 0.0

        return concentrations

    def __len__(self) -> int:
        """Return the number of nodes in this frame."""
        return len(self.nodes)

    def __str__(self) -> str:
        """Return a human-readable summary of the frame."""
        return f"Frame {self.frame_id} (num_nodes={len(self.nodes)}, num_clusters={len(self.clusters)})"

    def __repr__(self) -> str:
        """Return a detailed string representation of the frame."""
        return f"Frame {self.frame_id} (num_nodes={len(self.nodes)})\n(first node: {(self.nodes[0].symbol, self.nodes[0].position) if len(self.nodes) > 0 else ''}\n(lattice=\n{self.lattice})\n"

    def __del__(self) -> None:
        """Release references to internal data structures."""
        del self.nodes
        del self.clusters
        del self.lattice
        del self._data
        del self.connectivities
