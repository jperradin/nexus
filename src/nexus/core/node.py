import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List

from ..utils.geometry import wrap_position

@dataclass(slots=True, order=True)
class Node:
    """
    Dataclass representing a node in the simulation frame.

    A Node is the fundamental unit of the cluster analysis pipeline. It stores
    properties (symbol, position, mass) and participates in the union-find
    algorithm through its parent field. Nodes self-parent by default, meaning each
    node starts as the root of its own cluster.

    Attributes:
        symbol (str): Chemical symbol of the atom (e.g., "Si", "O").
        node_id (int): Unique identifier for the node.
        position (np.ndarray): 3D Cartesian coordinates of the atom.
        parent (Optional[Node]): Parent node in the union-find structure. Defaults to self,
            making each node its own root until merged via union().
        neighbors (List[Node]): List of neighboring nodes found by the neighbor searcher.
        cluster_id (Optional[int]): Identifier of the cluster this node belongs to.
            Defaults to node_id until assigned during clustering.
        distances (Optional[List[float]]): Distances to each neighbor, populated by the
            neighbor searcher.
        indices (Optional[List[int]]): Global indices of each neighbor in the frame's node list.
        mass (Optional[float]): Atomic, particle, or node mass in reduced units. Defaults to 1.0.
        coordination (Optional[int]): Coordination number (count of nearest neighbors).
            Defaults to 0.
        other (Optional[List[str]]): Additional per-atom attributes parsed from the trajectory file.
    """
    symbol: str
    node_id: int
    position: np.ndarray = field(compare=False, repr=False)
    parent: Optional['Node'] = field(default=None, compare=False, repr=False)
    neighbors: List['Node'] = field(default_factory=list, compare=False, repr=False)
    cluster_id: Optional[int] = field(default=None, compare=False, repr=False)
    distances: Optional[List[float]] = field(default=None, compare=False, repr=False)
    indices: Optional[List[int]] = field(default=None, compare=False, repr=False)
    mass: Optional[float] = field(default=None, compare=True, repr=True)
    coordination: Optional[int] = field(default=None, compare=True, repr=True)
    other: Optional[List[str]] = field(default=None, compare=False, repr=False)

    _next_id = 0

    def __post_init__(self):
        """
        Initialize default values after dataclass creation.

        Sets sensible defaults for optional fields that were not provided at construction
        time. Uses object.__setattr__ because the dataclass is configured with slots=True.
        The parent defaults to self to establish each node as its own union-find root.
        """

        if self.position is None:
            object.__setattr__(self, 'position', np.zeros(3))

        # Auto-increment node_id if not explicitly provided
        if self.node_id is None:
            object.__setattr__(self, 'node_id', Node._next_id)
            Node._next_id += 1

        if self.mass is None:
            object.__setattr__(self, 'mass', 1.0)

        if self.coordination is None:
            object.__setattr__(self, 'coordination', 0)

        if self.other is None:
            object.__setattr__(self, 'other', [])

        if self.neighbors is None:
            object.__setattr__(self, 'neighbors', [])

        # Self-parenting makes each node the root of its own cluster initially
        if self.parent is None:
            object.__setattr__(self, 'parent', self)

        if self.cluster_id is None:
            object.__setattr__(self, 'cluster_id', self.node_id)

    @staticmethod
    def wrap_position(position: np.ndarray, lattice: np.ndarray) -> np.ndarray:
        """
        Wrap a position back into the periodic simulation box.

        Delegates to the Numba-accelerated wrap_position utility function.

        Args:
            position (np.ndarray): 3D Cartesian coordinates to wrap.
            lattice (np.ndarray): 3x3 lattice matrix defining the simulation box.

        Returns:
            np.ndarray: The wrapped 3D coordinates inside the box.
        """
        return wrap_position(position, lattice)

    def add_neighbor(self, node: 'Node') -> None:
        """
        Append a node to this node's neighbor list.

        Args:
            node (Node): The neighboring node to add.
        """
        self.neighbors.append(node)

    def reset_parent(self) -> None:
        """
        Reset this node's parent to itself.

        Restores the node as the root of its own union-find set, effectively
        disconnecting it from any previously assigned cluster.
        """
        self.parent = self

    def set_coordination(self, coordination: int) -> None:
        """
        Set the coordination number for this node.

        Args:
            coordination (int): The number of nearest neighbors.
        """
        self.coordination = coordination

    def __str__(self) -> str:
        """Return a human-readable summary of the node."""
        return f"Node {self.node_id} ({self.symbol}) | Z = {self.coordination} | neighbors: {len(self.neighbors)} | position: {self.position}"

    def __repr__(self) -> str:
        """Return a string representation for debugging."""
        return self.__str__()
