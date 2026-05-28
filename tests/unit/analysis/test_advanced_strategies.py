"""Unit tests for the bond, coordination, and shared clustering strategies.

These run on a small Si-O bridged chain (Si-O-Si-O-Si-O-Si): 4 Si linked by 3
bridging O. With a Si-O cutoff of ~2.0 the bond criterion links all 4 Si into a
single Si-O-Si cluster; interior Si have O-coordination 2, terminal Si have 1.

get_connectivities() label generation (pairwise / mixing / alternating / default,
for both bond and distance criteria) is covered directly by mutating settings,
since it depends only on configuration, not geometry.
"""

import numpy as np
import pytest

from nexus.config.settings import (
    AnalysisSettings,
    ClusteringSettings,
    Cutoff,
    GeneralSettings,
    LatticeSettings,
    SettingsBuilder,
)
from nexus.analysis.strategies.bond_strategy import BondingStrategy
from nexus.analysis.strategies.coordination_strategy import CoordinationStrategy
from nexus.analysis.strategies.shared_strategy import SharedStrategy


# --------------------------------------------------------------------------- #
# Settings builders for the Si-O system
# --------------------------------------------------------------------------- #
def _general(path):
    return GeneralSettings(
        project_name="test",
        export_directory="test_export",
        file_location=path,
        range_of_frames=(0, -1),
        apply_pbc=True,
        verbose=False,
    )


def _cutoffs():
    return [
        Cutoff(type1="Si", type2="Si", distance=2.0),
        Cutoff(type1="Si", type2="O", distance=2.0),
        Cutoff(type1="O", type2="O", distance=2.0),
    ]


def bond_settings(path):
    clustering = ClusteringSettings(
        criterion="bond",
        node_types=["Si", "O"],
        node_masses=[28.0855, 15.9994],
        connectivity=["Si", "O", "Si"],
        cutoffs=_cutoffs(),
    )
    return _build(path, clustering)


def coordination_settings(path, **overrides):
    clustering = ClusteringSettings(
        criterion="bond",
        node_types=["Si", "O"],
        node_masses=[28.0855, 15.9994],
        connectivity=["Si", "O", "Si"],
        cutoffs=_cutoffs(),
        with_coordination_number=True,
        coordination_mode="O",
        coordination_range=[1, 2],
        with_connectivity_name="SiOSi",
    )
    for k, v in overrides.items():
        setattr(clustering, k, v)
    return _build(path, clustering)


def shared_settings(path):
    clustering = ClusteringSettings(
        criterion="bond",
        node_types=["Si", "O"],
        node_masses=[28.0855, 15.9994],
        connectivity=["Si", "O", "Si"],
        cutoffs=_cutoffs(),
        with_coordination_number=True,
        coordination_mode="O",
        coordination_range=[1, 2],
        with_number_of_shared=True,
        shared_mode="O",
        shared_threshold=1,
        shared_threshold_mode="exact",
        with_connectivity_name="corner",
    )
    return _build(path, clustering)


def _build(path, clustering):
    return (
        SettingsBuilder()
        .with_general(_general(path))
        .with_lattice(LatticeSettings(apply_custom_lattice=False))
        .with_clustering(clustering)
        .with_analysis(AnalysisSettings(with_all=True))
        .build()
    )


class TestBondStrategy:
    def test_links_all_si_into_one_cluster(self, sio2_chain_path, frame_factory):
        settings = bond_settings(sio2_chain_path)
        frame = frame_factory(sio2_chain_path, settings)
        strat = BondingStrategy(frame, settings)
        strat.find_neighbors()
        clusters = strat.build_clusters()
        assert len(clusters) == 1
        assert clusters[0].get_size() == 4
        assert clusters[0].get_connectivity() == "Si-O-Si"

    def test_connectivity_label(self, sio2_chain_path, frame_factory):
        settings = bond_settings(sio2_chain_path)
        frame = frame_factory(sio2_chain_path, settings)
        assert BondingStrategy(frame, settings).get_connectivities() == ["Si-O-Si"]

    def test_bad_connectivity_raises(self, sio2_chain_path, frame_factory):
        settings = bond_settings(sio2_chain_path)
        frame = frame_factory(sio2_chain_path, settings)
        settings.clustering.connectivity = ["Si", "Si"]  # only 2 -> invalid for bond
        with pytest.raises(ValueError):
            BondingStrategy(frame, settings).get_connectivities()


class TestCoordinationStrategy:
    def test_coordination_numbers_computed(self, sio2_chain_path, frame_factory):
        settings = coordination_settings(sio2_chain_path)
        frame = frame_factory(sio2_chain_path, settings)
        strat = CoordinationStrategy(frame, settings)
        strat.find_neighbors()
        # mode "O": interior Si bridge 2 oxygens, terminal Si bridge 1.
        si_coords = sorted(n.coordination for n in frame.get_nodes() if n.symbol == "Si")
        assert si_coords == [1, 1, 2, 2]

    def test_default_mode_builds_cluster(self, sio2_chain_path, frame_factory):
        settings = coordination_settings(sio2_chain_path)
        frame = frame_factory(sio2_chain_path, settings)
        strat = CoordinationStrategy(frame, settings)
        strat.find_neighbors()
        clusters = strat.build_clusters()
        assert len(clusters) == 1
        assert clusters[0].get_size() == 4
        assert clusters[0].get_connectivity() == "SiOSi"

    def test_calculate_coordination_modes(self, sio2_chain_path, frame_factory):
        # all_types counts every neighbor; interior Si has 2 O neighbors only here.
        settings = coordination_settings(sio2_chain_path, coordination_mode="all_types")
        frame = frame_factory(sio2_chain_path, settings)
        strat = CoordinationStrategy(frame, settings)
        strat.find_neighbors()
        si_coords = sorted(n.coordination for n in frame.get_nodes() if n.symbol == "Si")
        assert si_coords == [1, 1, 2, 2]


class TestSharedStrategy:
    def test_networking_nodes_are_interior_si(self, sio2_chain_path, frame_factory):
        settings = shared_settings(sio2_chain_path)
        frame = frame_factory(sio2_chain_path, settings)
        strat = SharedStrategy(frame, settings)
        strat.find_neighbors()
        # Only interior Si have 2 corner-sharing partners (terminal Si have 1).
        nodes = strat.get_networking_nodes(np.array([1, 2]))
        symbols = [n.symbol for n in nodes]
        assert symbols == ["Si", "Si"]

    def test_default_mode_builds_cluster(self, sio2_chain_path, frame_factory):
        settings = shared_settings(sio2_chain_path)
        frame = frame_factory(sio2_chain_path, settings)
        strat = SharedStrategy(frame, settings)
        strat.find_neighbors()
        clusters = strat.build_clusters()
        # Two interior Si share a bridging O -> one cluster of size 2.
        assert len(clusters) == 1
        assert clusters[0].get_size() == 2

    def test_distance_criterion_not_implemented(self, sio2_chain_path, frame_factory):
        settings = shared_settings(sio2_chain_path)
        frame = frame_factory(sio2_chain_path, settings)
        settings.clustering.criterion = "distance"
        strat = SharedStrategy(frame, settings)
        strat.find_neighbors()
        with pytest.raises(NotImplementedError):
            strat.build_clusters()


class TestGetConnectivitiesBranches:
    """Label generation across pairing modes, independent of geometry."""

    def _strat(self, sio2_chain_path, frame_factory, **overrides):
        settings = coordination_settings(sio2_chain_path, with_connectivity_name="", **overrides)
        frame = frame_factory(sio2_chain_path, settings)
        return CoordinationStrategy(frame, settings)

    def test_bond_pairwise(self, sio2_chain_path, frame_factory):
        strat = self._strat(sio2_chain_path, frame_factory, with_pairwise=True)
        assert strat.get_connectivities() == ["SiO_1-SiO_1", "SiO_2-SiO_2"]

    def test_bond_mixing(self, sio2_chain_path, frame_factory):
        strat = self._strat(sio2_chain_path, frame_factory, with_mixing=True)
        labels = strat.get_connectivities()
        # i in {1,2}, j>=i -> (1,1),(1,2),(2,2)
        assert labels == ["SiO_1-SiO_1", "SiO_1-SiO_2", "SiO_2-SiO_2"]

    def test_bond_alternating(self, sio2_chain_path, frame_factory):
        strat = self._strat(sio2_chain_path, frame_factory, with_alternating=True)
        assert len(strat.get_connectivities()) >= 2

    def test_distance_pairwise(self, sio2_chain_path, frame_factory):
        strat = self._strat(
            sio2_chain_path, frame_factory, criterion="distance", connectivity=["Si", "O"], with_pairwise=True
        )
        assert strat.get_connectivities() == ["Si_1-O_1", "Si_2-O_2"]

    def test_default_uses_connectivity_name(self, sio2_chain_path, frame_factory):
        settings = coordination_settings(sio2_chain_path)  # name "SiOSi", default mode
        frame = frame_factory(sio2_chain_path, settings)
        strat = CoordinationStrategy(frame, settings)
        assert strat.get_connectivities() == ["SiOSi"]
