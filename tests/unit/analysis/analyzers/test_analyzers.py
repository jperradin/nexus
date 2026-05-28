"""Unit tests for the analyzers in nexus.analysis.analyzers.

Each analyzer is fed a frame whose clusters were produced by the real
DistanceStrategy, so the per-frame -> finalize() pipeline is exercised end to
end. Two scenarios give known oracles:

  * two_clusters fixture -> 2 non-percolating clusters of size 3 each.
  * percolating fixture  -> 1 percolating cluster of size 27 (spans xyz).

Analyzers that exclude percolating clusters (average/spanning size, gyration,
correlation length, size distribution) therefore see the 2x3 clusters but
nothing in the percolating frame; percolation/order-parameter behave inversely.
"""

import numpy as np
import pytest

from nexus.analysis.analyzer_factory import AnalyzerFactory
from nexus.analysis.analyzers.average_cluster_size_analyzer import AverageClusterSizeAnalyzer
from nexus.analysis.analyzers.cluster_size_distribution_analyzer import ClusterSizeDistributionAnalyzer
from nexus.analysis.analyzers.concentration_analyzer import ConcentrationAnalyzer
from nexus.analysis.analyzers.correlation_length_analyzer import CorrelationLengthAnalyzer
from nexus.analysis.analyzers.gyration_radius_analyzer import GyrationRadiusAnalyzer
from nexus.analysis.analyzers.largest_cluster_size_analyzer import LargestClusterSizeAnalyzer
from nexus.analysis.analyzers.order_parameter_analyzer import OrderParameterAnalyzer
from nexus.analysis.analyzers.percolation_probability_analyzer import PercolationProbabilityAnalyzer
from nexus.analysis.analyzers.spanning_cluster_size_analyzer import SpanningClusterSizeAnalyzer
from nexus.analysis.strategies.distance_strategy import DistanceStrategy


def clustered_frame(path, settings_factory, frame_factory):
    """Cluster a fixture with DistanceStrategy and attach results to the frame."""
    settings = settings_factory(path)
    frame = frame_factory(path, settings)
    strat = DistanceStrategy(frame, settings)
    strat.find_neighbors()
    clusters = strat.build_clusters()
    conns = strat.get_connectivities()
    frame.set_clusters(clusters)
    frame.set_connectivities(conns)
    return frame, conns, settings


def run(analyzer, frame, conns):
    analyzer.analyze(frame, conns)
    return analyzer.finalize()


class TestLargestClusterSize:
    def test_percolating_lattice(self, percolating_path, settings_factory, frame_factory):
        frame, conns, settings = clustered_frame(percolating_path, settings_factory, frame_factory)
        res = run(LargestClusterSizeAnalyzer(settings), frame, conns)
        assert res["largest_cluster_size"]["A-A"] == 27

    def test_two_clusters(self, two_clusters_path, settings_factory, frame_factory):
        frame, conns, settings = clustered_frame(two_clusters_path, settings_factory, frame_factory)
        res = run(LargestClusterSizeAnalyzer(settings), frame, conns)
        assert res["largest_cluster_size"]["A-A"] == 3


class TestAverageClusterSize:
    def test_two_clusters_weight_average(self, two_clusters_path, settings_factory, frame_factory):
        # <S> = sum(s^2 n) / sum(s n) = (3^2 * 2) / (3 * 2) = 3.0
        frame, conns, settings = clustered_frame(two_clusters_path, settings_factory, frame_factory)
        res = run(AverageClusterSizeAnalyzer(settings), frame, conns)
        assert np.isclose(res["average_cluster_size"]["A-A"], 3.0)

    def test_percolating_excluded(self, percolating_path, settings_factory, frame_factory):
        # Only cluster is percolating -> excluded -> <S> = 0.
        frame, conns, settings = clustered_frame(percolating_path, settings_factory, frame_factory)
        res = run(AverageClusterSizeAnalyzer(settings), frame, conns)
        assert res["average_cluster_size"]["A-A"] == 0.0


class TestSpanningClusterSize:
    def test_two_clusters(self, two_clusters_path, settings_factory, frame_factory):
        frame, conns, settings = clustered_frame(two_clusters_path, settings_factory, frame_factory)
        res = run(SpanningClusterSizeAnalyzer(settings), frame, conns)
        assert res["spanning_cluster_size"]["A-A"] == 3

    def test_percolating_excluded(self, percolating_path, settings_factory, frame_factory):
        frame, conns, settings = clustered_frame(percolating_path, settings_factory, frame_factory)
        res = run(SpanningClusterSizeAnalyzer(settings), frame, conns)
        assert res["spanning_cluster_size"]["A-A"] == 0.0


class TestPercolationProbability:
    def test_percolating_is_one(self, percolating_path, settings_factory, frame_factory):
        frame, conns, settings = clustered_frame(percolating_path, settings_factory, frame_factory)
        res = run(PercolationProbabilityAnalyzer(settings), frame, conns)
        assert res["percolation_probabilities"]["A-A"] == 1.0

    def test_finite_is_zero(self, two_clusters_path, settings_factory, frame_factory):
        frame, conns, settings = clustered_frame(two_clusters_path, settings_factory, frame_factory)
        res = run(PercolationProbabilityAnalyzer(settings), frame, conns)
        assert res["percolation_probabilities"]["A-A"] == 0.0


class TestOrderParameter:
    def test_percolating_full(self, percolating_path, settings_factory, frame_factory):
        # All 27 networking nodes in the spanning cluster -> P_inf = 1.
        frame, conns, settings = clustered_frame(percolating_path, settings_factory, frame_factory)
        res = run(OrderParameterAnalyzer(settings), frame, conns)
        assert np.isclose(res["order_parameters"]["A-A"], 1.0)

    def test_finite_is_zero(self, two_clusters_path, settings_factory, frame_factory):
        frame, conns, settings = clustered_frame(two_clusters_path, settings_factory, frame_factory)
        res = run(OrderParameterAnalyzer(settings), frame, conns)
        assert res["order_parameters"]["A-A"] == 0.0


class TestConcentration:
    def test_all_networking_nodes_clustered(self, two_clusters_path, settings_factory, frame_factory):
        frame, conns, settings = clustered_frame(two_clusters_path, settings_factory, frame_factory)
        res = run(ConcentrationAnalyzer(settings), frame, conns)
        assert np.isclose(res["concentrations"]["A-A"], 1.0)


class TestGyrationRadius:
    def test_size_bin_present_and_positive(self, two_clusters_path, settings_factory, frame_factory):
        frame, conns, settings = clustered_frame(two_clusters_path, settings_factory, frame_factory)
        res = run(GyrationRadiusAnalyzer(settings), frame, conns)
        radii = res["gyration_radii"]["A-A"]
        # Clusters of size 3 -> a size-3 bin with a positive mean Rg.
        assert 3 in radii
        assert radii[3] > 0.0

    def test_percolating_empty(self, percolating_path, settings_factory, frame_factory):
        frame, conns, settings = clustered_frame(percolating_path, settings_factory, frame_factory)
        res = run(GyrationRadiusAnalyzer(settings), frame, conns)
        assert res["gyration_radii"]["A-A"] == {}


class TestCorrelationLength:
    def test_finite_positive(self, two_clusters_path, settings_factory, frame_factory):
        frame, conns, settings = clustered_frame(two_clusters_path, settings_factory, frame_factory)
        res = run(CorrelationLengthAnalyzer(settings), frame, conns)
        assert res["correlation_length"]["A-A"] > 0.0

    def test_percolating_zero(self, percolating_path, settings_factory, frame_factory):
        frame, conns, settings = clustered_frame(percolating_path, settings_factory, frame_factory)
        res = run(CorrelationLengthAnalyzer(settings), frame, conns)
        assert res["correlation_length"]["A-A"] == 0.0


class TestClusterSizeDistribution:
    def test_counts_two_clusters_of_size_three(self, two_clusters_path, settings_factory, frame_factory):
        frame, conns, settings = clustered_frame(two_clusters_path, settings_factory, frame_factory)
        res = run(ClusterSizeDistributionAnalyzer(settings), frame, conns)
        assert res["size_distribution"]["A-A"][3] == 2

    def test_percolating_empty(self, percolating_path, settings_factory, frame_factory):
        frame, conns, settings = clustered_frame(percolating_path, settings_factory, frame_factory)
        res = run(ClusterSizeDistributionAnalyzer(settings), frame, conns)
        assert res["size_distribution"]["A-A"] == {}


class TestAnalyzerFactory:
    def test_returns_registered_analyzer(self, default_settings):
        factory = AnalyzerFactory(default_settings, verbose=False)
        analyzer = factory.get_analyzer("LargestClusterSizeAnalyzer")
        assert isinstance(analyzer, LargestClusterSizeAnalyzer)

    def test_unknown_returns_none(self, default_settings):
        factory = AnalyzerFactory(default_settings, verbose=False)
        assert factory.get_analyzer("NopeAnalyzer") is None

    def test_all_nine_registered(self, default_settings):
        factory = AnalyzerFactory(default_settings, verbose=False)
        names = [
            "AverageClusterSizeAnalyzer",
            "ConcentrationAnalyzer",
            "LargestClusterSizeAnalyzer",
            "SpanningClusterSizeAnalyzer",
            "PercolationProbabilityAnalyzer",
            "OrderParameterAnalyzer",
            "ClusterSizeDistributionAnalyzer",
            "GyrationRadiusAnalyzer",
            "CorrelationLengthAnalyzer",
        ]
        assert all(factory.get_analyzer(n) is not None for n in names)
