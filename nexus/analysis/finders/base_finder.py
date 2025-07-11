from abc import ABC, abstractmethod
from typing import List

from ...core.frame import Frame
from ...core.node import Node
from ...config.settings import Settings


class BaseFinder(ABC):
    def __init__(self, frame: Frame, settings: Settings) -> None:
        self.frame = frame
        self._lattice = self.frame.lattice
        self._nodes = self.frame.nodes
        self._settings = settings

    @abstractmethod
    def find_neighbors(self) -> None:
        pass

    # common to all finders
    def find(self, node: Node) -> Node:
        if node.parent != node:
            node.parent = self.find(node.parent)
        return node.parent

    # common to all finders
    def union(self, node_1: Node, node_2: Node) -> None:
        root_1 = self.find(node_1)
        root_2 = self.find(node_2)
        
        if root_1 != root_2:
            root_2.parent = root_1

    @abstractmethod
    def find_clusters(self) -> None:
        pass

    @abstractmethod
    def get_connectivities(self) -> List[str]:
        pass

    def __str__(self) -> str:
        return f"{self.__class__.__name__}({self._settings})"

    def __repr__(self) -> str:
        return self.__str__()