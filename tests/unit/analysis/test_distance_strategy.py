"""Unit tests for nexus.analysis.strategies.distance_strategy.DistanceStrategy.

Exercises the full clustering pipeline (KD-tree neighbor search + union-find +
property calculation) on the hand-made fixtures whose expected cluster counts,
sizes and percolation flags are known by construction.

Note: build_clusters keeps only clusters of size > 1 (singletons are dropped).
"""

import numpy as np
import pytest

from nexus.analysis.strategies.distance_strategy import DistanceStrategy


def run_strategy(frame, settings):
    """Run the distance-clustering pipeline and return the cluster list."""
    strat = DistanceStrategy(frame, settings)
    strat.find_neighbors()
    return strat.build_clusters()


def sizes(clusters):
    return sorted(c.get_size() for c in clusters)


class TestConnectivityLabel:
    def test_two_element_label(self, two_atom_frame, default_settings):
        strat = DistanceStrategy(two_atom_frame, default_settings)
        assert strat.get_connectivities() == ["A-A"]

    def test_bad_connectivity_raises(self, two_atom_frame, settings_factory, two_atoms_path):
        settings = settings_factory(two_atoms_path)
        settings.clustering.connectivity = ["A", "A", "A"]  # invalid for distance
        strat = DistanceStrategy(two_atom_frame, settings)
        with pytest.raises(ValueError):
            strat.get_connectivities()


class TestClusterCounts:
    def test_two_atoms_form_one_cluster(self, two_atom_frame, default_settings):
        clusters = run_strategy(two_atom_frame, default_settings)
        assert len(clusters) == 1
        assert clusters[0].get_size() == 2

    def test_linear_chain_is_single_cluster(self, linear_chain_frame, linear_chain_path, settings_factory):
        clusters = run_strategy(linear_chain_frame, settings_factory(linear_chain_path))
        assert len(clusters) == 1
        assert clusters[0].get_size() == 5

    def test_two_disjoint_groups(self, two_cluster_frame, two_clusters_path, settings_factory):
        clusters = run_strategy(two_cluster_frame, settings_factory(two_clusters_path))
        assert len(clusters) == 2
        assert sizes(clusters) == [3, 3]


class TestPercolation:
    def test_full_lattice_percolates_all_axes(self, percolating_frame, percolating_path, settings_factory):
        clusters = run_strategy(percolating_frame, settings_factory(percolating_path))
        assert len(clusters) == 1
        cluster = clusters[0]
        assert cluster.get_size() == 27
        assert bool(cluster.is_percolating) is True
        assert set(cluster.percolation_probability) == set("xyz")


class TestPeriodicBoundary:
    def test_bond_only_across_boundary_with_pbc(self, pbc_pair_path, settings_factory, frame_factory):
        # x=0.25 and x=9.75, box 10: MIC separation 0.5 < cutoff -> bonded.
        settings = settings_factory(pbc_pair_path, apply_pbc=True)
        frame = frame_factory(pbc_pair_path, settings)
        clusters = run_strategy(frame, settings)
        assert len(clusters) == 1
        assert clusters[0].get_size() == 2

    def test_no_bond_without_pbc(self, pbc_pair_path, settings_factory, frame_factory):
        # Direct separation 9.5 > cutoff -> no bond -> both are singletons (dropped).
        settings = settings_factory(pbc_pair_path, apply_pbc=False)
        frame = frame_factory(pbc_pair_path, settings)
        clusters = run_strategy(frame, settings)
        assert len(clusters) == 0


class TestClusterProperties:
    def test_concentration_and_order_parameter_set(self, percolating_frame, percolating_path, settings_factory):
        clusters = run_strategy(percolating_frame, settings_factory(percolating_path))
        cluster = clusters[0]
        # Single spanning cluster: every networking node belongs to it.
        assert np.isclose(cluster.concentration, 1.0)
        # Percolates in 3 dims -> order parameter non-zero on all axes.
        assert all(p > 0.0 for p in cluster.order_parameter)
