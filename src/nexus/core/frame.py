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
    Reprensation of a frame of a trajectory

    Attributes:
    -----------
    frame_id : int
        Id of the frame
    nodes : List[Node]
        List of nodes in the frame
    lattice : np.ndarray
        Lattice of the frame
    clusters : Optional[List[Cluster]]
        List of clusters in the frame
    _data : Dict[str, np.ndarray]
        Internal data structure for node data (symbol, position)
    """

    frame_id: int
    nodes: List[Node]
    lattice: np.ndarray
    _data: Dict[str, np.ndarray]
    _settings: Settings
    clusters: Optional[List[Cluster]] = None
    connectivities: Optional[List[str]] = None

    def __post_init__(self):
        """Initialisation after object creation"""
        if not isinstance(self.nodes, list):
            raise TypeError("nodes must be a list of Nodes")
        if self.lattice is not None and not isinstance(self.lattice, np.ndarray):
            raise TypeError("lattice must be a numpy array")

    def initialize_nodes(self) -> None:
        """Initialize the list of nodes in the frame, filters out not selected species"""
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
        """Set the lattice of the frame"""
        if lattice.shape != (3, 3):
            raise ValueError("lattice must be a 3x3 numpy array")

        try:
            np.linalg.inv(lattice)
        except np.linalg.LinAlgError:
            raise ValueError("lattice must be a non-singular matrix")

        self.lattice = lattice

    def get_lattice(self) -> Optional[np.ndarray]:
        """Get the lattice of the frame"""
        return self.lattice

    def get_unique_elements(self) -> List[str]:
        """Get the unique elements in the frame"""
        return np.unique([node.symbol for node in self.nodes])

    def get_node_by_id(self, node_id: int) -> Optional[Node]:
        """Get an node by its id"""
        for node in self.nodes:
            if node.node_id == node_id:
                return node
        return None

    def get_positions(self) -> np.ndarray:
        """Get the positions of all nodes in the frame"""
        return np.array([node.position for node in self.nodes])

    def get_positions_by_element(self) -> Dict[str, np.ndarray]:
        """Get the positions of all nodes in the frame grouped by element"""
        return {
            node.symbol: np.array(
                [node.position for node in self.nodes if node.symbol == node.symbol]
            )
            for node in self.nodes
        }

    def get_wrapped_positions(self) -> np.ndarray:
        """Get the wrapped positions of all nodes in the frame"""
        return wrap_positions(self.get_positions(), self.lattice)

    def get_wrapped_positions_by_element(self) -> Dict[str, np.ndarray]:
        """Get the wrapped positions of all nodes in the frame grouped by element"""
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
        """Get the clusters of the frame"""
        return self.clusters

    def get_nodes(self) -> List[Node]:
        """Get the nodes of the frame"""
        return self.nodes

    def get_networking_nodes(self) -> int:
        total_sizes = [c.get_size() for c in self.clusters]
        return np.sum(total_sizes)

    def get_connectivities(self) -> List[str]:
        """Get the connectivities of the frame"""
        return self.connectivities

    def set_connectivities(self, connectivities: List[str]) -> None:
        """Set the connectivities of the frame"""
        self.connectivities = connectivities

    def add_cluster(self, cluster: Cluster) -> None:
        """Add a cluster to the frame"""
        if self.clusters is None:
            self.clusters = []
        self.clusters.append(cluster)

    def set_clusters(self, clusters: List[Cluster]) -> None:
        """Set the clusters of the frame"""
        for cluster in clusters:
            cluster.frame_id = self.frame_id
            cluster.set_lattice(self.lattice)
        self.clusters = clusters

    def get_concentration(self) -> Dict[str, float]:
        """Get the concentrations of each cluster connectivity in the frame"""
        concentrations = {}
        for connectivity in self.connectivities:
            for cluster in self.clusters:
                if cluster.get_connectivity() == connectivity:
                    concentrations[connectivity] = cluster.total_nodes / len(self.nodes)
                    break

        for connectivity in self.connectivities:
            if connectivity not in concentrations:
                concentrations[connectivity] = 0.0

        return concentrations

    def __len__(self) -> int:
        """Get the number of nodes in the frame"""
        return len(self.nodes)

    def __str__(self) -> str:
        return f"Frame {self.frame_id} (num_nodes={len(self.nodes)}, num_clusters={len(self.clusters)})"

    def __repr__(self) -> str:
        return f"Frame {self.frame_id} (num_nodes={len(self.nodes)})\n(first node: {(self.nodes[0].symbol, self.nodes[0].position) if len(self.nodes) > 0 else ''}\n(lattice=\n{self.lattice})\n"

    def __del__(self) -> None:
        del self.nodes
        del self.clusters
        del self.lattice
        del self._data
        del self.connectivities

