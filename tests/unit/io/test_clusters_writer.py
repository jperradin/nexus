"""Unit tests for nexus.io.writer.clusters_writer.ClustersWriter.

Distance clusters (single type A) cover the print-mode dispatch and XYZ output;
a bond Si-O-Si cluster with unwrapped-cluster printing enabled covers the
decoration-node and .bonds output paths (bonds are only written for the bond
criterion).
"""

import numpy as np

from nexus.analysis.strategies.bond_strategy import BondingStrategy
from nexus.analysis.strategies.distance_strategy import DistanceStrategy
from nexus.config.settings import (
    AnalysisSettings,
    ClusteringSettings,
    Cutoff,
    GeneralSettings,
    LatticeSettings,
    SettingsBuilder,
)
from nexus.io.writer.clusters_writer import ClustersWriter


def distance_clusters(path, settings_factory, frame_factory, tmp_path, print_mode):
    settings = settings_factory(path, export_directory=str(tmp_path))
    settings.clustering.print_mode = print_mode
    frame = frame_factory(path, settings)
    strat = DistanceStrategy(frame, settings)
    strat.find_neighbors()
    clusters = strat.build_clusters()
    frame.set_clusters(clusters)  # assigns frame_id + lattice to each cluster
    return settings, frame.get_clusters()


def bond_clusters(path, frame_factory, tmp_path):
    clustering = ClusteringSettings(
        criterion="bond",
        node_types=["Si", "O"],
        node_masses=[28.0855, 15.9994],
        connectivity=["Si", "O", "Si"],
        cutoffs=[
            Cutoff("Si", "Si", 2.0),
            Cutoff("Si", "O", 2.0),
            Cutoff("O", "O", 2.0),
        ],
        with_printed_unwrapped_clusters=True,
        print_mode="connectivity",
    )
    settings = (
        SettingsBuilder()
        .with_general(
            GeneralSettings(
                project_name="test",
                export_directory=str(tmp_path),
                file_location=path,
                range_of_frames=(0, -1),
                apply_pbc=True,
                verbose=False,
            )
        )
        .with_lattice(LatticeSettings(apply_custom_lattice=False))
        .with_clustering(clustering)
        .with_analysis(AnalysisSettings(with_all=True))
        .build()
    )
    frame = frame_factory(path, settings)
    strat = BondingStrategy(frame, settings)
    strat.find_neighbors()
    clusters = strat.build_clusters()
    frame.set_clusters(clusters)
    return settings, frame.get_clusters()


class TestPrintModeDispatch:
    def test_none_writes_nothing(self, percolating_path, settings_factory, frame_factory, tmp_path):
        settings, clusters = distance_clusters(percolating_path, settings_factory, frame_factory, tmp_path, "none")
        w = ClustersWriter(settings)
        w.set_clusters(clusters)
        w.write()
        assert not (tmp_path / "unwrapped_clusters").exists()

    def test_all_mode_writes_single_xyz(self, percolating_path, settings_factory, frame_factory, tmp_path):
        settings, clusters = distance_clusters(percolating_path, settings_factory, frame_factory, tmp_path, "all")
        w = ClustersWriter(settings)
        w.set_clusters(clusters)
        w.write()
        xyz = tmp_path / "unwrapped_clusters" / "all_unwrapped_clusters-frame_0.xyz"
        assert xyz.exists()
        assert int(xyz.read_text().splitlines()[0]) == 27

    def test_connectivity_mode_subdir(self, percolating_path, settings_factory, frame_factory, tmp_path):
        settings, clusters = distance_clusters(percolating_path, settings_factory, frame_factory, tmp_path, "connectivity")
        w = ClustersWriter(settings)
        w.set_clusters(clusters)
        w.write()
        assert (tmp_path / "unwrapped_clusters" / "A-A" / "A-A_unwrapped_clusters-frame_0.xyz").exists()

    def test_individual_mode_per_cluster_file(self, two_clusters_path, settings_factory, frame_factory, tmp_path):
        settings, clusters = distance_clusters(two_clusters_path, settings_factory, frame_factory, tmp_path, "individual")
        w = ClustersWriter(settings)
        w.set_clusters(clusters)
        w.write()
        files = list((tmp_path / "unwrapped_clusters" / "A-A").glob("cluster-frame_0-id_*.xyz"))
        assert len(files) == 2

    def test_sorted_by_size_descending(self, two_clusters_path, settings_factory, frame_factory, tmp_path):
        settings, clusters = distance_clusters(two_clusters_path, settings_factory, frame_factory, tmp_path, "none")
        w = ClustersWriter(settings)
        w.set_clusters(clusters)
        sizes = [c.size for c in w._clusters]
        assert sizes == sorted(sizes, reverse=True)


class TestBondOutput:
    def test_writes_xyz_and_bonds(self, sio2_chain_path, frame_factory, tmp_path):
        settings, clusters = bond_clusters(sio2_chain_path, frame_factory, tmp_path)
        w = ClustersWriter(settings)
        w.set_clusters(clusters)
        w.write()
        base = tmp_path / "unwrapped_clusters" / "Si-O-Si"
        xyz = base / "Si-O-Si_unwrapped_clusters-frame_0.xyz"
        bonds = base / "Si-O-Si_unwrapped_clusters-frame_0.bonds"
        assert xyz.exists() and bonds.exists()
        # Bonds file lists Si-Si linkages for the bond criterion.
        assert "Si(" in bonds.read_text()
        # XYZ includes bridging O decoration atoms (4 Si + 3 O = 7).
        assert int(xyz.read_text().splitlines()[0]) == 7
