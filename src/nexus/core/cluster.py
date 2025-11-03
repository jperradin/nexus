from typing import List, Tuple, Set, Dict
import numpy as np
import os
from scipy.sparse.csgraph import breadth_first_order
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

    #### Old method which is not accurate enough
    # def calculate_percolation_probability(self) -> None:
    #     """
    #     Pre-determine percolation probability based on gyration radius and lattice dimensions.
    #     The cluster percolates if it spans all three dimensions.
    #     """
    #     if self.size <= 1: return
    #     min_coords = np.min(self.unwrapped_positions, axis=0)
    #     max_coords = np.max(self.unwrapped_positions, axis=0)
    #     span = max_coords - min_coords
    #     percolate_x = span[0] > self.lattice[0, 0]
    #     percolate_y = span[1] > self.lattice[1, 1]
    #     percolate_z = span[2] > self.lattice[2, 2]
    #     self.percolation_probability = ""
    #     if percolate_x: self.percolation_probability += 'x'
    #     if percolate_y: self.percolation_probability += 'y'
    #     if percolate_z: self.percolation_probability += 'z'
    #     self.is_percolating = 'x' in self.percolation_probability and 'y' in self.percolation_probability and 'z' in self.percolation_probability

    def calculate_percolation_probability(self) -> None:
        """
        Robust percolation detection using period vector algorithm.
        Detects true percolation by checking if cluster connects to itself
        across periodic boundaries through actual bond connectivity.

        Method suggested by Livraghi et al. (2021) https://pubs.acs.org/action/showCitFormats?doi=10.1021/acs.jctc.1c00423&ref=pdf

        """
        if self.size <= 1:
            self.is_percolating = False
            self.percolation_probability = ""
            return
        
        # Detect period vectors during BFS traversal
        period_vectors = self._detect_period_vectors()
        
        # Calculate percolation dimension from period vectors
        percolation_dim = self._calculate_period_dimension(period_vectors)
        
        # Determine which directions percolate
        self.percolation_probability = self._get_percolation_directions(period_vectors)
        
        # True 3D percolation requires spanning all three dimensions
        self.is_percolating = percolation_dim == 3

    def _detect_period_vectors(self) -> List[np.ndarray]:
        """
        Detect period vectors by traversing the cluster and finding when
        we encounter the same node through different periodic paths.
        """
        if not self.nodes:
            return []
        
        root_node = self.nodes[0].parent
        queue = [root_node]
        
        # Track unwrapped positions during traversal
        visited_positions = {root_node.node_id: root_node.position.copy()}
        visited_set = {root_node.node_id}
        period_vectors = []
        
        if self.settings.clustering.criterion == 'distance':
            while queue:
                current_node = queue.pop(0)
                current_unwrapped = visited_positions[current_node.node_id]
                
                for neighbor in current_node.neighbors:
                    # Only consider neighbors in this cluster
                    if neighbor.cluster_id != self.root_id:
                        continue
                    
                    # Calculate unwrapped position of neighbor
                    relative_vec = self._unwrap_vector(neighbor.position - current_node.position)
                    neighbor_unwrapped = current_unwrapped + relative_vec
                    
                    if neighbor.node_id in visited_set:
                        # Already visited - check for period vector
                        period = neighbor_unwrapped - visited_positions[neighbor.node_id]
                        
                        # If period is non-zero, we found a periodic connection
                        if np.linalg.norm(period) > 1e-6:
                            period_vectors.append(period)
                    else:
                        # First time visiting this neighbor
                        visited_positions[neighbor.node_id] = neighbor_unwrapped
                        visited_set.add(neighbor.node_id)
                        queue.append(neighbor)
        elif self.settings.clustering.criterion == 'bond':
            bridge_symbol = self.settings.clustering.connectivity[1]
            while queue:
                current_node = queue.pop(0)
                current_unwrapped = visited_positions[current_node.node_id]
                
                for neighbor in current_node.neighbors:
                    if neighbor.symbol == bridge_symbol:
                        for neighbor2 in neighbor.neighbors:
                            # Only consider neighbors in this cluster
                            if neighbor2.cluster_id != self.root_id:
                                continue
                            
                            # Calculate unwrapped position of neighbor
                            relative_vec = self._unwrap_vector(neighbor2.position - current_node.position)
                            neighbor_unwrapped = current_unwrapped + relative_vec
                            
                            if neighbor2.node_id in visited_set:
                                # Already visited - check for period vector
                                period = neighbor_unwrapped - visited_positions[neighbor2.node_id]
                                
                                # If period is non-zero, we found a periodic connection
                                if np.linalg.norm(period) > 1e-6:
                                    period_vectors.append(period)
                            else:
                                # First time visiting this neighbor
                                visited_positions[neighbor2.node_id] = neighbor_unwrapped
                                visited_set.add(neighbor2.node_id)
                                queue.append(neighbor2)
        
        return period_vectors

    def _calculate_period_dimension(self, period_vectors: List[np.ndarray]) -> int:
        """
        Calculate algebraic dimension of period vector set using SVD.
        Returns 0, 1, 2, or 3 indicating percolation dimensionality.
        """
        if not period_vectors:
            return 0
        
        # Stack period vectors into matrix
        period_matrix = np.array(period_vectors)
        
        # Use SVD to find rank (number of linearly independent vectors)
        _, singular_values, _ = np.linalg.svd(period_matrix)
        
        # Count non-zero singular values (with numerical tolerance)
        tolerance = 1e-6 * singular_values[0] if len(singular_values) > 0 else 1e-10
        rank = np.sum(singular_values > tolerance)
        
        return min(rank, 3)  # Cap at 3 for 3D systems

    def _get_percolation_directions(self, period_vectors: List[np.ndarray]) -> str:
        """
        Determine which Cartesian directions (x, y, z) the cluster percolates in
        based on linearly independent period vectors.
        """
        if not period_vectors:
            return ""
        
        # Find linearly independent period vectors
        independent_periods = self._get_independent_periods(period_vectors)
        
        # Check which directions each independent period spans
        directions = {'x': False, 'y': False, 'z': False}
        lattice_diag = np.diag(self.lattice)
        
        for period in independent_periods:
            # Normalize by lattice vectors to get fractional coordinates
            fractional = period / lattice_diag
            
            # Check if period has significant component in each direction
            if abs(fractional[0]) > 0.5:
                directions['x'] = True
            if abs(fractional[1]) > 0.5:
                directions['y'] = True
            if abs(fractional[2]) > 0.5:
                directions['z'] = True
        
        result = ""
        if directions['x']: result += 'x'
        if directions['y']: result += 'y'
        if directions['z']: result += 'z'
        
        return result

    def _get_independent_periods(self, period_vectors: List[np.ndarray]) -> List[np.ndarray]:
        """
        Extract linearly independent period vectors using Gram-Schmidt-like process.
        """
        if not period_vectors:
            return []
        
        independent = []
        tolerance = 1e-8
        
        for period in period_vectors:
            if np.linalg.norm(period) < tolerance:
                continue
            
            if not independent:
                independent.append(period)
            else:
                # Check if adding this vector increases the span
                test_matrix = np.vstack(independent + [period])
                rank_before = np.linalg.matrix_rank(np.vstack(independent), tol=tolerance)
                rank_after = np.linalg.matrix_rank(test_matrix, tol=tolerance)
                
                if rank_after > rank_before:
                    independent.append(period)
                
                # Stop when we have 3 independent vectors (max for 3D)
                if len(independent) == 3:
                    break
        return independent

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

        if self.settings.clustering.criterion == 'distance':
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
        elif self.settings.clustering.criterion == 'bond':
            bridge_symbol = self.settings.clustering.connectivity[1]
            while queue:
                current_node = queue.pop(0)
                for neighbor in current_node.neighbors:
                    if neighbor.symbol == bridge_symbol:
                        for neighbor2 in neighbor.neighbors:
                            if neighbor2.cluster_id == self.root_id and neighbor2.node_id not in visited_nodes:
                                relative_position = self._unwrap_vector(neighbor2.position - current_node.position)
                                dict_positions[neighbor2.node_id] = dict_positions[current_node.node_id] + relative_position
                                link = tuple(sorted((current_node.node_id, neighbor2.node_id)))
                                self._linkage_set.add(link)
                                visited_nodes.add(neighbor2.node_id)
                                queue.append(neighbor2)
                                pbar.update(1)
        
        pbar.close()
        self.set_indices_and_positions(dict_positions)
        self.linkages = sorted(list(self._linkage_set))
        
        # --- Find and unwrap decorating nodes ---
        if self.settings.clustering.criterion == 'bond' and self.settings.clustering.with_printed_unwrapped_clusters:
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
