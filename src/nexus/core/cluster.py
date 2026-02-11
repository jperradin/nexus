from typing import List, Tuple, Set, Dict
import numpy as np
from tqdm import tqdm
import shutil

# Internal imports
from .node import Node
from ..config.settings import Settings
from ..utils.geometry import calculate_gyration_radius, wrap_position


class Cluster:
    """
    Represents a group of connected nodes forming a single cluster.

    A cluster is produced by a clustering strategy and consumed by analyzers. It holds
    the set of member nodes and provides methods to compute physical properties such as
    unwrapped positions, center of mass, gyration radius, percolation probability, and
    order parameter.

    Attributes:
        nodes (List[Node]): Member nodes belonging to this cluster.
        connectivity (str): Connectivity label identifying the cluster type (e.g., "Si-Si").
        root_id (int): Identifier of the root node used as the cluster ID.
        size (int): Number of member nodes in the cluster.
        settings (Settings): Configuration settings used during analysis.
        lattice (np.ndarray): 3x3 lattice matrix of the simulation cell.
        center_of_mass (np.ndarray): Wrapped center-of-mass position.
        symbols (list): Chemical symbols of member nodes in insertion order.
        indices (list): Node IDs of member nodes in insertion order.
        unwrapped_positions (np.ndarray): Cartesian positions unwrapped across periodic
            boundaries, shape (N, 3).
        percolation_probability (str): String of percolating directions (e.g., "xyz", "xz").
        gyration_radius (float): Radius of gyration computed from unwrapped positions.
        order_parameter (list): Per-dimension order parameter [P_x, P_y, P_z].
        total_nodes (int): Total number of nodes in the frame, set externally by analyzers.
        concentration (float): Fraction of total nodes belonging to this cluster.
        is_percolating (bool): True if the cluster percolates in all three dimensions.
        is_spanning (bool): True if the cluster spans the simulation cell.
        linkages (List[Tuple[int, int]]): Sorted list of unique node-pair bonds in the cluster.
        period_vectors (List[np.ndarray]) : List of periodic vectors of the cluster, empty if non-periodic.
        decoration_atoms (Dict[int, Dict]): Bridging nodes not part of the cluster but
            connected to member nodes, keyed by node ID.
    """

    def __init__(
        self,
        connectivity: str,
        root_id: int,
        size: int,
        settings: Settings,
        lattice: np.ndarray,
    ) -> None:
        """
        Initialize a cluster with its connectivity label, root node, size, and lattice.

        Args:
            connectivity (str): Connectivity label identifying the cluster type.
            root_id (int): Identifier of the root node, also used as the cluster ID.
            size (int): Number of member nodes in the cluster.
            settings (Settings): Configuration settings for the analysis.
            lattice (np.ndarray): 3x3 lattice matrix of the simulation cell.
        """
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
        self.percolation_probability: str = ""
        self.gyration_radius: float = 0.0
        self.order_parameter: list = [0.0] * 3
        self.total_nodes: int = 0
        self.concentration: float = 0.0
        self.is_percolating: bool = False
        self.is_spanning: bool = False

        self.linkages: List[Tuple[int, int]] = []
        self._linkage_set: Set[Tuple[int, int]] = set()
        self.period_vectors: List[np.ndarray] = []

        # New attribute to store decorating nodes (e.g., bridging oxygens)
        self.decoration_atoms: Dict[int, Dict] = {}

    def add_node(self, node: Node) -> None:
        """
        Add a node to this cluster and assign it the cluster ID.

        Args:
            node (Node): The node to add. Its ``cluster_id`` is set to this cluster's
                ``root_id``.
        """
        node.cluster_id = self.root_id
        self.nodes.append(node)

    def set_lattice(self, lattice: np.ndarray) -> None:
        """
        Update the lattice matrix and recompute its inverse.

        Args:
            lattice (np.ndarray): New 3x3 lattice matrix of the simulation cell.
        """
        self.lattice = lattice
        self._inv_lattice = np.linalg.inv(lattice)

    def get_nodes(self) -> List[Node]:
        """
        Return the list of member nodes.

        Returns:
            List[Node]: Nodes belonging to this cluster.
        """
        return self.nodes

    def get_connectivity(self) -> str:
        """
        Return the connectivity label of this cluster.

        Returns:
            str: Connectivity label (e.g., "Si-Si" or "Si-O-Si").
        """
        return self.connectivity

    def get_size(self) -> int:
        """
        Return the number of member nodes in this cluster.

        Returns:
            int: Cluster size.
        """
        return self.size

    def set_indices_and_positions(self, positions_dict) -> None:
        """
        Populate symbols, indices, and unwrapped positions from a positions dictionary.

        Iterates over member nodes and matches them against the provided dictionary to
        build the ``symbols``, ``indices``, and ``unwrapped_positions`` arrays.

        Args:
            positions_dict (dict): Mapping of node IDs to their unwrapped Cartesian
                positions.
        """
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
        """
        Compute the wrapped center of mass and translate unwrapped positions accordingly.

        Calculates the mean of the unwrapped positions, wraps it back into the simulation
        cell, and shifts all unwrapped positions (including decoration nodes) by the
        resulting translation vector so they are centered around the wrapped center of mass.
        """
        if self.unwrapped_positions.size > 0:
            unwrapped_com = np.mean(self.unwrapped_positions, axis=0)
            wrapped_com = wrap_position(unwrapped_com, self.lattice)
            translation_vector = wrapped_com - unwrapped_com
            self.unwrapped_positions += translation_vector
            # Also translate decoration nodes
            for atom_id in self.decoration_atoms:
                self.decoration_atoms[atom_id]["position"] += translation_vector
            self.center_of_mass = wrapped_com
        else:
            self.center_of_mass = np.zeros(3)

    def calculate_gyration_radius(self) -> None:
        """
        Compute the radius of gyration from unwrapped positions.

        Delegates the heavy computation to the Numba-accelerated utility function.
        Computes the center of mass first if it has not been calculated yet.
        """
        if self.size <= 1 or self.unwrapped_positions.size == 0:
            self.gyration_radius = 0.0
            return
        if not np.any(self.center_of_mass):
            self.calculate_center_of_mass()
        self.gyration_radius = calculate_gyration_radius(
            self.unwrapped_positions, self.center_of_mass
        )

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
        Detect percolation using the period vector algorithm.

        Determines whether the cluster connects to itself across periodic boundaries
        through actual bond connectivity, rather than relying on spatial extent alone.
        The algorithm proceeds in three steps: (1) detect period vectors via BFS traversal,
        (2) compute the percolation dimensionality from linearly independent period vectors,
        and (3) map period vectors to Cartesian directions.

        Reference:
            Livraghi et al. (2021)
            https://pubs.acs.org/action/showCitFormats?doi=10.1021/acs.jctc.1c00423&ref=pdf
        """
        if self.size <= 1:
            self.is_percolating = False
            self.percolation_probability = ""
            return

        # Detect period vectors during BFS traversal
        # period_vectors = self._detect_period_vectors()

        # Calculate percolation dimension from period vectors
        # percolation_dim = self._calculate_period_dimension(period_vectors)
        percolation_dim = self._calculate_period_dimension()

        # Determine which directions percolate
        self.percolation_probability = self._get_percolation_directions()

        # True 3D percolation requires spanning all three dimensions
        self.is_percolating = percolation_dim == 3

    def _calculate_period_dimension(self) -> int:
        """
        Calculate the algebraic dimension of a set of period vectors using SVD.

        Stacks the period vectors into a matrix and computes its numerical rank via
        singular value decomposition. The rank equals the number of linearly independent
        periodic directions, capped at 3.

        Args:
            period_vectors (List[np.ndarray]): Period vectors detected during BFS traversal.

        Returns:
            int: Percolation dimensionality (0, 1, 2, or 3).
        """
        if not self.period_vectors:
            return 0

        # Stack period vectors into matrix
        period_matrix = np.array(self.period_vectors)

        # Use SVD to find rank (number of linearly independent vectors)
        _, singular_values, _ = np.linalg.svd(period_matrix)

        # Count non-zero singular values (with numerical tolerance)
        tolerance = 1e-6 * singular_values[0] if len(singular_values) > 0 else 1e-10
        rank = np.sum(singular_values > tolerance)

        return min(rank, 3)  # Cap at 3 for 3D systems

    def _get_percolation_directions(self) -> str:
        """
        Determine which Cartesian directions the cluster percolates in.

        Extracts linearly independent period vectors, converts them to fractional
        coordinates, and checks whether each direction has a significant component
        (greater than half a lattice vector).

        Args:
            period_vectors (List[np.ndarray]): Period vectors detected during BFS traversal.

        Returns:
            str: Concatenation of percolating direction labels (e.g., "x", "xy", "xyz").
        """
        if not self.period_vectors:
            return ""

        # Find linearly independent period vectors
        independent_periods = self._get_independent_periods(self.period_vectors)
        if not independent_periods:
            return ""

        # Check which directions each independent period spans
        directions = {"x": False, "y": False, "z": False}
        lattice_diag = np.diag(self.lattice)

        for period in independent_periods:
            # Normalize by lattice vectors to get fractional coordinates
            fractional = period / lattice_diag

            # Check if period has significant component in each direction
            if abs(fractional[0]) > 0.5:
                directions["x"] = True
            if abs(fractional[1]) > 0.5:
                directions["y"] = True
            if abs(fractional[2]) > 0.5:
                directions["z"] = True

        result = ""
        if directions["x"]:
            result += "x"
        if directions["y"]:
            result += "y"
        if directions["z"]:
            result += "z"

        return result

    def _get_independent_periods(
        self, period_vectors: List[np.ndarray]
    ) -> List[np.ndarray]:
        """
        Extract linearly independent period vectors using an incremental rank test.

        Iterates over period vectors and retains those that increase the rank of the
        accumulated matrix, up to a maximum of 3 independent vectors.

        Args:
            period_vectors (List[np.ndarray]): Candidate period vectors to filter.

        Returns:
            List[np.ndarray]: Subset of linearly independent period vectors (at most 3).
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
                rank_before = np.linalg.matrix_rank(
                    np.vstack(independent), tol=tolerance
                )
                rank_after = np.linalg.matrix_rank(test_matrix, tol=tolerance)

                if rank_after > rank_before:
                    independent.append(period)

                # Stop when we have 3 independent vectors (max for 3D)
                if len(independent) == 3:
                    break
        return independent

    def calculate_order_parameter(self) -> None:
        """
        Compute the per-dimension order parameter.

        Calculates P_inf as the ratio of cluster size to total nodes and assigns it to
        each percolating dimension. Non-percolating dimensions receive a value of zero.
        """
        if self.size <= 1 or self.total_nodes == 0:
            return
        p_inf = self.size / self.total_nodes
        if len(self.percolation_probability) == 1:
            self.order_parameter = [p_inf, 0.0, 0.0]
        elif len(self.percolation_probability) == 2:
            self.order_parameter = [p_inf, p_inf, 0.0]
        elif len(self.percolation_probability) == 3:
            self.order_parameter = [p_inf, p_inf, p_inf]

    def calculate_concentration(self) -> None:
        """
        Compute the concentration of this cluster as the ratio of its size to total nodes.
        """
        if self.total_nodes > 0:
            self.concentration = self.size / self.total_nodes

    def _unwrap_vector(self, vector: np.ndarray) -> np.ndarray:
        """
        Apply minimum-image convention to a displacement vector.

        Converts the vector to fractional coordinates, rounds to the nearest image,
        and converts back to Cartesian coordinates.

        Args:
            vector (np.ndarray): Cartesian displacement vector.

        Returns:
            np.ndarray: Minimum-image displacement vector in Cartesian coordinates.
        """
        fractional_vector = np.dot(vector, self._inv_lattice)
        fractional_vector -= np.round(fractional_vector)
        return np.dot(fractional_vector, self.lattice)

    def calculate_unwrapped_positions(self) -> None:
        """
        Unwrap node positions across periodic boundaries using BFS traversal.

        Starting from the root node, traverses the cluster graph and reconstructs
        continuous (unwrapped) positions by applying the minimum-image convention to each
        neighbor displacement. Tracks unique node-pair linkages during traversal. For the
        bond criterion, also unwraps decoration nodes (bridging nodes) connected to member
        nodes.

        """
        # TODO : stores period_vectors on the fly for percolation probability calculation
        if self.size <= 1:
            return

        root_node = self.nodes[0].parent
        queue = [root_node]
        visited_positions = {root_node.node_id: root_node.position}
        visited_nodes = {root_node.node_id}

        # periodic vector for percolation probability calculation
        period_vectors = []

        progress_bar_kwargs = {
            "disable": not self.settings.verbose,
            "leave": False,
            "ncols": shutil.get_terminal_size().columns,
            "colour": "magenta",
        }
        pbar_desc = f"Unwrapping cluster {self.root_id} ({self.connectivity})"
        pbar = tqdm(total=self.size, desc=pbar_desc, **progress_bar_kwargs)
        pbar.update(1)

        # Distance criterion: direct neighbor-to-neighbor traversal
        if self.settings.clustering.criterion == "distance":
            while queue:
                current_node = queue.pop(0)
                current_unwrapped = visited_positions[current_node.node_id]

                for neighbor in current_node.neighbors:
                    # Only consider neighbors in this cluster
                    if neighbor.cluster_id != self.root_id:
                        continue

                    # Calculate unwrapped position of neighbor
                    relative_vec = self._unwrap_vector(
                        neighbor.position - current_node.position
                    )

                    neighbor_unwrapped = current_unwrapped + relative_vec

                    # Store period vector if already visited
                    if neighbor.node_id in visited_nodes:
                        period = (
                            neighbor_unwrapped - visited_positions[neighbor.node_id]
                        )

                        # If period is non-zero -> found a periodic connection
                        if np.linalg.norm(period) > 1e-6:
                            period_vectors.append(period)
                    else:
                        # First time visiting this neighbor
                        visited_positions[neighbor.node_id] = neighbor_unwrapped
                        link = tuple(sorted((current_node.node_id, neighbor.node_id)))
                        self._linkage_set.add(link)
                        visited_nodes.add(neighbor.node_id)
                        queue.append(neighbor)
                        pbar.update(1)
            
        # Bond criterion: traverse through bridging nodes (e.g., Si -> O -> Si)
        elif self.settings.clustering.criterion == "bond":
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
                            relative_vec = self._unwrap_vector(
                                neighbor2.position - current_node.position
                            )
                            neighbor_unwrapped = current_unwrapped + relative_vec

                            if neighbor2.node_id in visited_nodes:
                                # Already visited - check for period vector
                                period = (
                                    neighbor_unwrapped - visited_positions[neighbor2.node_id]
                                )

                                # If period is non-zero -> found a periodic connection
                                if np.linalg.norm(period) > 1e-6:
                                    period_vectors.append(period)
                            else:
                                # First time visiting this neighbor
                                visited_positions[neighbor2.node_id] = neighbor_unwrapped
                                link = tuple(sorted((current_node.node_id, neighbor2.node_id)))
                                self._linkage_set.add(link)
                                visited_nodes.add(neighbor2.node_id)
                                queue.append(neighbor2)
                                pbar.update(1)

        pbar.close()
        self.set_indices_and_positions(visited_positions)
        self.linkages = sorted(list(self._linkage_set))
        self.period_vectors = period_vectors

        # --- Find and unwrap decorating nodes ---
        if (
            self.settings.clustering.criterion == "bond"
            and self.settings.clustering.with_printed_unwrapped_clusters
        ):
            bridge_symbol = self.settings.clustering.connectivity[1]
            # Create a map of node ID to original Node object for quick lookup
            cluster_nodes_map = {node.node_id: node for node in self.nodes}

            for node_id, unwrapped_pos in visited_positions.items():
                original_node = cluster_nodes_map.get(node_id)
                if not original_node:
                    continue

                for neighbor in original_node.neighbors:
                    if (
                        neighbor.symbol == bridge_symbol
                        and neighbor.node_id not in self.decoration_atoms
                    ):
                        relative_pos = self._unwrap_vector(
                            neighbor.position - original_node.position
                        )
                        unwrapped_bridge_pos = unwrapped_pos + relative_pos
                        self.decoration_atoms[neighbor.node_id] = {
                            "symbol": neighbor.symbol,
                            "position": unwrapped_bridge_pos,
                            "coordination": neighbor.coordination,
                        }

    def __str__(self) -> str:
        """Return a human-readable summary of the cluster."""
        list_id = [str(i.node_id) for i in self.nodes]
        if len(list_id) > 20:
            list_id = list_id[:20] + ["..."]
        return f"{self.root_id} {self.connectivity} {self.size} {self.is_percolating} {', '.join(list_id)}"

    def __repr__(self) -> str:
        """Return the string representation of the cluster."""
        return self.__str__()
