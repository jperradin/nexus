import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List

from ..utils.geometry import wrap_position

@dataclass(slots=True, order=True)
class Node:
    """ Reprensation of a node 
    
    Attributes:
    -----------
    symbol : str
        Symbol of the node
    node_id : int
        Id of the node (unique identifier, autoincremented)
    position : np.ndarray
        Position of the node
    parent : Optional['Node']
        Parent of the node (Optional)
    neighbors : List['Node']
        List of neighbors of the node (Optional)
    distances : Optional[List[float]]
        Distances of the neighbors of the node (Optional)
    indices : Optional[List[int]]
        Indices of the neighbors of the node (Optional)
    mass : float
        Mass of the node (Optional)
    coordination : Optional[int]
        Coordination number of the node (Optional)
    other : Optional[List[str]]
        Other attributes of the node (Optional)
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
        """ Initialisation after object creation """

        if self.position is None:
            object.__setattr__(self, 'position', np.zeros(3))
        
        if self.node_id is None:
            object.__setattr__(self, 'node_id', Node._next_id)
            Node._next_id += 1
        
        if self.mass is None:
            object.__setattr__(self, 'mass', 0.0)      

        if self.coordination is None:
            object.__setattr__(self, 'coordination', 0)      

        if self.other is None:
            object.__setattr__(self, 'other', [])

        if self.neighbors is None:
            object.__setattr__(self, 'neighbors', [])

        if self.parent is None:
            object.__setattr__(self, 'parent', self)

        if self.cluster_id is None:
            object.__setattr__(self, 'cluster_id', self.node_id)    

    @staticmethod
    def wrap_position(position: np.ndarray, lattice: np.ndarray) -> np.ndarray:
        """ Wrap position in a periodic box defined by the lattice """
        return wrap_position(position, lattice)

    def add_neighbor(self, node: 'Node') -> None:
        """ Add a node as a neighbor """
        self.neighbors.append(node)

    def reset_parent(self) -> None:
        self.parent = self

    def set_coordination(self, coordination: int) -> None:
        self.coordination = coordination

    def __str__(self) -> str:
        return f"Node {self.node_id} ({self.symbol}) | coordination: {self.coordination} | neighbors: {len(self.neighbors)} | position: {self.position}"

    def __repr__(self) -> str:
        return self.__str__()
