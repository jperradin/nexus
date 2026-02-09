import numpy as np
from scipy.spatial import cKDTree
from tqdm import tqdm
from typing import List
import shutil

from ....core.node import Node
from ....core.frame import Frame
from ....config.settings import Settings
from ....utils.geometry import cartesian_to_fractional


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
        """Build the KD-tree and assign neighbors to every node in the frame."""
        positions = self.frame.get_wrapped_positions()

        # Build the k-d tree, handling periodic boundary conditions
        if self.settings.apply_pbc:
            positions_frac = cartesian_to_fractional(positions, self._lattice)
            kdtree = cKDTree(positions_frac, boxsize=[1, 1, 1])
            query_positions = positions_frac
            # Estimate fractional cutoff. This is an approximation but is only used for broad-phase search.
            # The exact distance check will perform the precise filtering.
            search_radius = (
                self._max_cutoff / np.linalg.norm(self._lattice, axis=0).max()
            )
        else:
            kdtree = cKDTree(positions)
            query_positions = positions
            search_radius = self._max_cutoff

        progress_bar_kwargs = {
            "disable": not self.settings.verbose,
            "leave": False,
            "ncols": shutil.get_terminal_size().columns,
            "colour": "green",
        }

        progress_bar = tqdm(
            range(len(self._nodes)),
            desc="Fetching nearest neighbors ...",
            **progress_bar_kwargs,
        )

        for i in progress_bar:
            node = self._nodes[i]

            # Find candidate neighbors within the max cutoff radius
            indices = kdtree.query_ball_point(query_positions[i], search_radius)

            # Refine neighbors with exact distance checks
            self._filter_and_assign_neighbors(node, indices)

    def _filter_and_assign_neighbors(
        self, node: Node, candidate_indices: List[int]
    ) -> None:
        """
        Refine candidates with exact per-pair cutoff checks and assign results.

        Args:
            node (Node): The node whose neighbors are being assigned.
            candidate_indices (List[int]): Indices from the broad-phase KD-tree
                query.
        """
        new_neighbors = []
        new_distances = []

        node_pos = node.position

        for neighbor_idx in candidate_indices:
            neighbor = self._nodes[neighbor_idx]

            # Skip self-interaction
            if node.node_id == neighbor.node_id:
                continue

            # Check exact cutoff distance for this pair of node types
            rcut = self.settings.clustering.get_cutoff(node.symbol, neighbor.symbol)
            if rcut is None:
                continue

            # Calculate distance (PBC or direct)
            if self.settings.apply_pbc:
                # You should implement calculate_pbc_distance in geometry.py
                # For now, we assume it exists
                from ....utils.geometry import calculate_pbc_distance

                dist = calculate_pbc_distance(
                    node_pos, neighbor.position, self._lattice
                )
            else:
                dist = np.linalg.norm(node_pos - neighbor.position)

            if dist <= rcut:
                new_neighbors.append(neighbor)
                new_distances.append(dist)

        node.neighbors = new_neighbors
        node.distances = new_distances
        node.indices = [n.node_id for n in new_neighbors]
