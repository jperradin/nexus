import numpy as np
from scipy.spatial import cKDTree
from tqdm import tqdm
from typing import List
import shutil

from ....core.node import Node
from ....core.frame import Frame
from ....config.settings import Settings
from ....utils.geometry import (
    cartesian_to_fractional,
    filter_neighbors_pbc_batch,
    filter_neighbors_direct_batch,
)


class NeighborSearcher:
    """
    KD-tree based neighbor finder for all nodes in a frame.

    Builds a ``cKDTree`` over node positions (with optional PBC support),
    queries it for candidates within the largest cutoff, then refines
    results with exact per-pair distance checks.

    Attributes:
        frame (Frame): The frame containing the nodes to process.
        settings (Settings): Configuration settings with cutoffs and PBC flag.
        _nodes (List[Node]): Direct reference to the nodes in the frame.
        _lattice (np.ndarray): Lattice matrix from the frame.
        _max_cutoff (float): Largest cutoff distance across all pair types.
    """

    def __init__(self, frame: Frame, settings: Settings):
        """
        Initialize the neighbor searcher.

        Args:
            frame (Frame): The frame containing the nodes to process.
            settings (Settings): Configuration settings.
        """
        self.frame: Frame = frame
        self.settings: Settings = settings
        self._nodes: List[Node] = frame.nodes
        self._lattice: np.ndarray = frame.lattice
        self._max_cutoff: float = max(
            c.distance for c in self.settings.clustering.cutoffs
        )

    def execute(self) -> None:
        positions = self.frame.get_wrapped_positions()
        N = positions.shape[0]

        lattice = np.ascontiguousarray(self._lattice, dtype=np.float64)
        inv_lattice = np.ascontiguousarray(np.linalg.inv(lattice), dtype=np.float64)

        # KD-tree (broad-phase)
        if self.settings.apply_pbc:
            positions_frac = cartesian_to_fractional(positions, lattice)
            kdtree = cKDTree(positions_frac, boxsize=[1.0, 1.0, 1.0])
            query_positions = positions_frac
            search_radius = (
                self._max_cutoff / np.linalg.norm(lattice, axis=0).max()
            )
        else:
            kdtree = cKDTree(positions)
            query_positions = positions
            search_radius = self._max_cutoff

        # Per-pair cutoff lookup table (squared)
        sym_to_id, rcut2_table = self._build_cutoff_table()

        # Symbol id per node (np.int64)
        sym_ids = np.empty(N, dtype=np.int64)
        unknown_id = -1
        for i, n in enumerate(self._nodes):
            sym_ids[i] = sym_to_id.get(n.symbol, unknown_id)

        # Cartesian positions array (separate from query_positions which may be fractional)
        cart_positions = np.ascontiguousarray(positions, dtype=np.float64)

        # Bulk broad-phase query (parallel where supported)
        try:
            all_candidates = kdtree.query_ball_point(
                query_positions, search_radius, workers=-1
            )
        except TypeError:
            all_candidates = kdtree.query_ball_point(query_positions, search_radius)

        progress_bar_kwargs = {
            "disable": not self.settings.verbose,
            "leave": False,
            "ncols": shutil.get_terminal_size().columns,
            "colour": "green",
        }
        progress_bar = tqdm(
            range(N),
            desc="Fetching nearest neighbors ...",
            **progress_bar_kwargs,
        )

        apply_pbc = self.settings.apply_pbc

        for i in progress_bar:
            cand_idx = all_candidates[i]
            if not cand_idx:
                self._assign_empty(self._nodes[i])
                continue

            # Drop self-interaction
            cand_arr = np.fromiter(
                (j for j in cand_idx if j != i), dtype=np.int64
            )
            if cand_arr.size == 0:
                self._assign_empty(self._nodes[i])
                continue

            # Build per-candidate squared cutoffs via table lookup
            si = sym_ids[i]
            if si < 0:
                self._assign_empty(self._nodes[i])
                continue
            cand_syms = sym_ids[cand_arr]
            rcut2 = rcut2_table[si, cand_syms]  # negative entries = invalid pair

            cand_pos = cart_positions[cand_arr]

            if apply_pbc:
                keep, dists = filter_neighbors_pbc_batch(
                    cart_positions[i], cand_pos, lattice, inv_lattice, rcut2
                )
            else:
                keep, dists = filter_neighbors_direct_batch(
                    cart_positions[i], cand_pos, rcut2
                )

            kept_global = cand_arr[keep]
            kept_dists = dists[keep]

            node = self._nodes[i]
            node.neighbors = [self._nodes[j] for j in kept_global]
            node.distances = kept_dists.tolist()
            node.indices = [self._nodes[j].node_id for j in kept_global]

    @staticmethod
    def _assign_empty(node: Node) -> None:
        node.neighbors = []
        node.distances = []
        node.indices = []

    def _build_cutoff_table(self):
        """
        Build (symbol -> int_id) map and a 2D squared-cutoff table.

        Returns:
            sym_to_id (Dict[str, int])
            rcut2_table (np.ndarray): (S, S) float64. Entries are squared
                cutoff distances; ``-1.0`` flags an invalid pair.
        """
        symbols: List[str] = []
        seen = set()
        # Stable ordering: settings.clustering.node_types first, then any extra
        # symbols seen in cutoffs.
        for s in self.settings.clustering.node_types:
            if s not in seen:
                seen.add(s); symbols.append(s)
        for c in self.settings.clustering.cutoffs:
            for s in (c.type1, c.type2):
                if s not in seen:
                    seen.add(s); symbols.append(s)

        sym_to_id: Dict[str, int] = {s: i for i, s in enumerate(symbols)}
        S = len(symbols)
        table = np.full((S, S), -1.0, dtype=np.float64)
        for c in self.settings.clustering.cutoffs:
            i = sym_to_id[c.type1]
            j = sym_to_id[c.type2]
            d2 = float(c.distance) * float(c.distance)
            table[i, j] = d2
            table[j, i] = d2  # symmetric
        return sym_to_id, table

    ##### Old implementation, slower.
    # def execute(self) -> None:
    #     """Build the KD-tree and assign neighbors to every node in the frame."""
    #     positions = self.frame.get_wrapped_positions()
    #
    #     # Build the k-d tree, handling periodic boundary conditions
    #     if self.settings.apply_pbc:
    #         positions_frac = cartesian_to_fractional(positions, self._lattice)
    #         kdtree = cKDTree(positions_frac, boxsize=[1, 1, 1])
    #         query_positions = positions_frac
    #         # Estimate fractional cutoff. This is an approximation but is only used for broad-phase search.
    #         # The exact distance check will perform the precise filtering.
    #         search_radius = (
    #             self._max_cutoff / np.linalg.norm(self._lattice, axis=0).max()
    #         )
    #     else:
    #         kdtree = cKDTree(positions)
    #         query_positions = positions
    #         search_radius = self._max_cutoff
    #
    #     progress_bar_kwargs = {
    #         "disable": not self.settings.verbose,
    #         "leave": False,
    #         "ncols": shutil.get_terminal_size().columns,
    #         "colour": "green",
    #     }
    #
    #     progress_bar = tqdm(
    #         range(len(self._nodes)),
    #         desc="Fetching nearest neighbors ...",
    #         **progress_bar_kwargs,
    #     )
    #
    #     for i in progress_bar:
    #         node = self._nodes[i]
    #
    #         # Find candidate neighbors within the max cutoff radius
    #         indices = kdtree.query_ball_point(query_positions[i], search_radius)
    #
    #         # Refine neighbors with exact distance checks
    #         self._filter_and_assign_neighbors(node, indices)
    #
    # def _filter_and_assign_neighbors(
    #     self, node: Node, candidate_indices: List[int]
    # ) -> None:
    #     """
    #     Refine candidates with exact per-pair cutoff checks and assign results.
    #
    #     Args:
    #         node (Node): The node whose neighbors are being assigned.
    #         candidate_indices (List[int]): Indices from the broad-phase KD-tree
    #             query.
    #     """
    #     new_neighbors = []
    #     new_distances = []
    #
    #     node_pos = node.position
    #
    #     for neighbor_idx in candidate_indices:
    #         neighbor = self._nodes[neighbor_idx]
    #
    #         # Skip self-interaction
    #         if node.node_id == neighbor.node_id:
    #             continue
    #
    #         # Check exact cutoff distance for this pair of node types
    #         rcut = self.settings.clustering.get_cutoff(node.symbol, neighbor.symbol)
    #         if rcut is None:
    #             continue
    #
    #         # Calculate distance (PBC or direct)
    #         if self.settings.apply_pbc:
    #             # You should implement calculate_pbc_distance in geometry.py
    #             # For now, we assume it exists
    #             from ....utils.geometry import calculate_pbc_distance
    #
    #             dist = calculate_pbc_distance(
    #                 node_pos, neighbor.position, self._lattice
    #             )
    #         else:
    #             dist = np.linalg.norm(node_pos - neighbor.position)
    #
    #         if dist <= rcut:
    #             new_neighbors.append(neighbor)
    #             new_distances.append(dist)
    #
    #     node.neighbors = new_neighbors
    #     node.distances = new_distances
    #     node.indices = [n.node_id for n in new_neighbors]
