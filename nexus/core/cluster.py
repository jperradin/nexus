from typing import List
import numpy as np
import os
from tqdm import tqdm

# Internal imports
from .node import Node
from ..config.settings import Settings

class Cluster:
    def __init__(self, connectivity: str, root_id: int, size: int, settings: Settings, lattice: np.ndarray) -> None:
        self.nodes: List[Node] = []
        self.connectivity: str = connectivity
        self.root_id: int = root_id
        self.size: int = size
        self.settings: Settings = settings
        self.lattice: np.ndarray = lattice

        self.center_of_mass: list = []
        self.symbols: list = []
        self.indices: list = []
        self.unwrapped_positions: list = []
        self.percolation_probability: str = ''
        self.gyration_radius: float = 0.0
        self.order_parameter: list = [0.0] * 3
        self.total_nodes: int = 0
        self.concentration: float = 0.0
        self.is_percolating: bool = False
        
    def add_node(self, node: Node) -> None:
        node.cluster_id = self.root_id
        self.nodes.append(node)

    def set_lattice(self, lattice: np.ndarray) -> None:
        self.lattice = lattice

    def get_nodes(self) -> List[Node]:
        return self.nodes

    def get_connectivity(self) -> str:
        return self.connectivity

    def get_size(self) -> int:
        return self.size

    def set_indices_and_positions(self, positions_dict) -> None:
        for node_id, position in positions_dict.items():
            for node in self.nodes:
                if node.node_id == node_id:
                    self.symbols.append(node.symbol)
                    break
            self.indices.append(node_id)
            self.unwrapped_positions.append([position[0], position[1], position[2]])
        self.unwrapped_positions = np.array(self.unwrapped_positions)

    def calculate_center_of_mass(self) -> None:
        self.center_of_mass = np.mean(self.unwrapped_positions, axis=0)
        
    def calculate_gyration_radius(self) -> None:
        self.gyration_radius = 0.0
        if self.size <= 1:
            return
        
        for i in range(self.unwrapped_positions.shape[0]):
            squared_rij = np.linalg.norm(self.unwrapped_positions[i, :] - self.unwrapped_positions[:, :], axis=1)** 2
            self.gyration_radius += np.sum(squared_rij)
            
        # Normalize the sum by 0.5 s²
        self.gyration_radius = np.sqrt((0.5 / (self.size**2)) * self.gyration_radius) 

    def calculate_percolation_probability(self) -> None:
        percolate_x = False
        percolate_y = False
        percolate_z = False

        if self.size <= 1:
            return

        for i in range(self.unwrapped_positions.shape[0]):
            dx = np.abs(self.unwrapped_positions[i, 0] - self.unwrapped_positions[:, 0])
            dy = np.abs(self.unwrapped_positions[i, 1] - self.unwrapped_positions[:, 1])
            dz = np.abs(self.unwrapped_positions[i, 2] - self.unwrapped_positions[:, 2])
            
            dx = np.max(dx)
            dy = np.max(dy)
            dz = np.max(dz)
            
            if dx > self.lattice[0, 0]:
                percolate_x = True
            if dy > self.lattice[1, 1]:
                percolate_y = True
            if dz > self.lattice[2, 2]:
                percolate_z = True
        
        if percolate_x:
            self.percolation_probability += 'x'
        if percolate_y:
            self.percolation_probability += 'y'
        if percolate_z:
            self.percolation_probability += 'z'
        
        # True if the cluster is percolating in all directions
        self.is_percolating = 'x' in self.percolation_probability and 'y' in self.percolation_probability and 'z' in self.percolation_probability
        # self.is_percolating = 'x' in self.percolation_probability or 'y' in self.percolation_probability or 'z' in self.percolation_probability

    def calculate_order_parameter(self) -> None:
        if self.size <= 1:
            return
        elif len(self.percolation_probability) == 1:
            self.order_parameter[0] = self.size / self.total_nodes
            self.order_parameter[1] = 0.0
            self.order_parameter[2] = 0.0
        elif len(self.percolation_probability) == 2:
            self.order_parameter[0] = self.size / self.total_nodes
            self.order_parameter[1] = self.size / self.total_nodes
            self.order_parameter[2] = 0.0
        elif len(self.percolation_probability) == 3:
            self.order_parameter[0] = self.size / self.total_nodes
            self.order_parameter[1] = self.size / self.total_nodes
            self.order_parameter[2] = self.size / self.total_nodes

    def calculate_concentration(self) -> None:
        self.concentration = self.size / self.total_nodes

    def calculate_unwrapped_positions(self) -> None:
        stack = [self.nodes[0].parent]

        if self.size <= 1:
            return

        dict_positions = {stack[0].node_id: self.nodes[0].position}

        progress_bar_kwargs = {
            "disable": not self.settings.verbose,
            "leave": False,
            "ncols": os.get_terminal_size().columns,
            "colour": "magenta"
        }

        criteria = self.settings.clustering.criteria
        chain = self.settings.clustering.connectivity

        if criteria == 'bond':
            node_1 = chain[0]
            bridge = chain[1]
            node_2 = chain[2]
            while tqdm(stack, desc=f"Unwrapping clusters {self.root_id} {self.connectivity} ...", **progress_bar_kwargs):
                current_node = stack.pop()
                if current_node.symbol == node_1:
                    for fn in current_node.neighbors:
                        if fn.symbol == bridge:
                            for sn in fn.neighbors:
                                if (
                                    sn.symbol == node_2 
                                    and sn.node_id not in dict_positions
                                    and sn.cluster_id == self.root_id
                                ):
                                    # Compute relative position from the current atom to its second_neighbour
                                    relative_position = self.unwrap_position(
                                        sn.position - current_node.position
                                    )

                                    # Accumulate relative position to get unwrapped position
                                    dict_positions[sn.node_id] = (
                                        dict_positions[current_node.node_id] + relative_position
                                    )

                                    # Add second_neighbour to the stack
                                    stack.append(sn)
        elif criteria == 'distance':
            node_1 = chain[0]
            node_2 = chain[1]
            while tqdm(stack, desc=f"Unwrapping clusters {self.root_id} {self.connectivity} ...", **progress_bar_kwargs):
                current_node = stack.pop()
                if current_node.symbol == node_1:
                    for fn in current_node.neighbors:
                        if (
                            fn.symbol == node_2
                            and fn.node_id not in dict_positions
                            and fn.cluster_id == self.root_id
                        ):
                            # Compute relative position from the current atom to its first neighbour
                            relative_position = self.unwrap_position(
                                fn.position - current_node.position
                            )
                            
                            # Accumulate relative position to get unwrapped position
                            dict_positions[fn.node_id] = (
                                dict_positions[current_node.node_id] + relative_position
                            )
                            
                            # Add first neighbour to the stack
                            stack.append(fn)
        
        self.set_indices_and_positions(dict_positions)

    def unwrap_position(self, vector):
        """
        Unwraps position considering periodic boundary conditions.
        """
        unwrapped_position = []

        for i in range(3):
            delta = vector[i] - round(vector[i] / self.lattice[i, i]) * self.lattice[i, i]
            unwrapped_position.append(delta)
        return tuple(unwrapped_position)

    def __str__(self) -> str:
        list_id = [str(i.node_id) for i in self.nodes]
        if len(list_id) > 20:
            list_id = list_id[:20] + ['...']
        list_id = ', '.join(list_id)
            
        return f"{self.root_id} {self.connectivity} {self.size} {self.is_percolating} {list_id}"

    def __repr__(self) -> str:
        return self.__str__()
        