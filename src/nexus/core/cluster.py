from typing import List, Tuple, Set, Dict
import numpy as np
import os
from tqdm import tqdm

# Internal imports
from .node import Node
from ..config.settings import Settings
from ..utils.geometry import calculate_gyration_radius, wrap_position

class Cluster:
    def __init__(self, connectivity: str, root_id: int, size: int, settings: Settings, lattice: np.ndarray) -> None:
        self.nodes: List[Node] = []
        self.connectivity: str = connectivity
        self.root_id: int = root_id
        self.size: int = size
        self.settings: Settings = settings
        self.lattice: np.ndarray = lattice
        self._inv_lattice: np.ndarray = np.linalg.inv(lattice)
        self._all_connectivities: Set[str] = set()

        self.center_of_mass: np.ndarray = np.zeros(3)
        self.symbols: list = []
        self.indices: list = []
        self.unwrapped_positions: np.ndarray = np.array([])
        self.percolation_probability: str = ''
        self.gyration_radius: float = 0.0
        self.order_parameter: list = [0.0] * 3
        self.total_nodes: int = 0
        self.concentration: float = 0.0
        self.is_percolating: bool = False
        self.is_spanning: bool = False
        
        self.linkages: List[Tuple[int, int]] = []
        self._linkage_set: Set[Tuple[int, int]] = set()
        
        # New attribute to store decorating nodes (e.g., bridging oxygens)
        self.decoration_atoms: Dict[int, Dict] = {}

    def add_node(self, node: Node) -> None:
        node.cluster_id = self.root_id
        self.nodes.append(node)

    def set_lattice(self, lattice: np.ndarray) -> None:
        self.lattice = lattice
        self._inv_lattice = np.linalg.inv(lattice)

    def get_nodes(self) -> List[Node]:
        return self.nodes

    def get_connectivity(self) -> str:
        return self.connectivity

    def get_size(self) -> int:
        return self.size

    def set_indices_and_positions(self, positions_dict) -> None:
        unwrapped_pos_list = []
        for node in self.nodes:
            for node_id, position in positions_dict.items():
                if node.node_id == node_id:
                    self.symbols.append(node.symbol)
                    self.indices.append(node.node_id)
                    unwrapped_pos_list.append(position)
                    break
        self.unwrapped_positions = np.array(unwrapped_pos_list)

    def calculate_center_of_mass(self) -> None:
        if self.unwrapped_positions.size > 0:
            unwrapped_com = np.mean(self.unwrapped_positions, axis=0)
            wrapped_com = wrap_position(unwrapped_com, self.lattice)
            translation_vector = wrapped_com - unwrapped_com
            self.unwrapped_positions += translation_vector
            # Also translate decoration nodes
            for atom_id in self.decoration_atoms:
                self.decoration_atoms[atom_id]['position'] += translation_vector
            self.center_of_mass = wrapped_com
        else:
            self.center_of_mass = np.zeros(3)
        
    def calculate_gyration_radius(self) -> None:
        if self.size <= 1 or self.unwrapped_positions.size == 0:
            self.gyration_radius = 0.0
            return
        if not np.any(self.center_of_mass):
            self.calculate_center_of_mass()
        self.gyration_radius = calculate_gyration_radius(self.unwrapped_positions, self.center_of_mass)

    def calculate_percolation_probability(self) -> None:
        """
        Pre-determine percolation probability based on gyration radius and lattice dimensions.
        The cluster percolates if it spans all three dimensions.
        """
        if self.size <= 1: return
        min_coords = np.min(self.unwrapped_positions, axis=0)
        max_coords = np.max(self.unwrapped_positions, axis=0)
        span = max_coords - min_coords
        percolate_x = span[0] > self.lattice[0, 0]
        percolate_y = span[1] > self.lattice[1, 1]
        percolate_z = span[2] > self.lattice[2, 2]
        self.percolation_probability = ""
        if percolate_x: self.percolation_probability += 'x'
        if percolate_y: self.percolation_probability += 'y'
        if percolate_z: self.percolation_probability += 'z'
        self.is_percolating = 'x' in self.percolation_probability and 'y' in self.percolation_probability and 'z' in self.percolation_probability

    def verify_percolation(self) -> None:
        """ 
        Verify if cluster truly percolates in given directory by checking
        connectivity between opposite boundary slices through periodic images. 

        This function is called after calculate_percolation_probability.
        """
        if not self.is_percolating:
            return

        axis_mapping = {'x': 0, 'y': 1, 'z': 2}
        slices = {axis: [] for axis in axis_mapping}
        # Determine atoms in each slice        
        slice_thickness = self.settings.clustering.get_max_cutoff() * 1.5
        


    def calculate_order_parameter(self) -> None:
        if self.size <= 1 or self.total_nodes == 0: return
        p_inf = self.size / self.total_nodes
        if len(self.percolation_probability) == 1: self.order_parameter = [p_inf, 0.0, 0.0]
        elif len(self.percolation_probability) == 2: self.order_parameter = [p_inf, p_inf, 0.0]
        elif len(self.percolation_probability) == 3: self.order_parameter = [p_inf, p_inf, p_inf]

    def calculate_concentration(self) -> None:
        if self.total_nodes > 0: self.concentration = self.size / self.total_nodes

    def _unwrap_vector(self, vector: np.ndarray) -> np.ndarray:
        fractional_vector = np.dot(vector, self._inv_lattice)
        fractional_vector -= np.round(fractional_vector)
        return np.dot(fractional_vector, self.lattice)

    def calculate_unwrapped_positions(self) -> None:
        if self.size <= 1: return

        root_node = self.nodes[0].parent
        queue = [root_node]
        dict_positions = {root_node.node_id: root_node.position}
        visited_nodes = {root_node.node_id}

        progress_bar_kwargs = {"disable": not self.settings.verbose, "leave": False, "ncols": os.get_terminal_size().columns, "colour": "magenta"}
        pbar_desc = f"Unwrapping cluster {self.root_id} ({self.connectivity})"
        pbar = tqdm(total=self.size, desc=pbar_desc, **progress_bar_kwargs)
        pbar.update(1)

        while queue:
            current_node = queue.pop(0)
            for neighbor in current_node.neighbors:
                if neighbor.cluster_id == self.root_id and neighbor.node_id not in visited_nodes:
                    relative_position = self._unwrap_vector(neighbor.position - current_node.position)
                    dict_positions[neighbor.node_id] = dict_positions[current_node.node_id] + relative_position
                    link = tuple(sorted((current_node.node_id, neighbor.node_id)))
                    self._linkage_set.add(link)
                    visited_nodes.add(neighbor.node_id)
                    queue.append(neighbor)
                    pbar.update(1)
        
        pbar.close()
        self.set_indices_and_positions(dict_positions)
        self.linkages = sorted(list(self._linkage_set))
        
        # --- New logic to find and unwrap decorating nodes ---
        if self.settings.clustering.criterion == 'bond':
            bridge_symbol = self.settings.clustering.connectivity[1]
            # Create a map of node ID to original Node object for quick lookup
            cluster_nodes_map = {node.node_id: node for node in self.nodes}

            for node_id, unwrapped_pos in dict_positions.items():
                original_node = cluster_nodes_map.get(node_id)
                if not original_node: continue

                for neighbor in original_node.neighbors:
                    if neighbor.symbol == bridge_symbol and neighbor.node_id not in self.decoration_atoms:
                        relative_pos = self._unwrap_vector(neighbor.position - original_node.position)
                        unwrapped_bridge_pos = unwrapped_pos + relative_pos
                        self.decoration_atoms[neighbor.node_id] = {
                            'symbol': neighbor.symbol,
                            'position': unwrapped_bridge_pos,
                            'coordination': neighbor.coordination
                        }

    def __str__(self) -> str:
        list_id = [str(i.node_id) for i in self.nodes]
        if len(list_id) > 20: list_id = list_id[:20] + ['...']
        return f"{self.root_id} {self.connectivity} {self.size} {self.is_percolating} {', '.join(list_id)}"

    def __repr__(self) -> str:
        return self.__str__()
