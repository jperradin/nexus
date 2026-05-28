"""Shared fixtures for the nexus-cat test suite.

Provides:
  - Path helpers to the hand-made trajectory fixtures in ``tests/data/``.
  - ``Settings`` factories for the distance-clustering pipeline (single type "A").
  - ``Frame`` factories that parse a fixture file the same way the real reader does.
  - Lightweight ``Node`` and ``Cluster`` builders for pure unit tests.

All geometric fixtures use known coordinates so expected values (distances,
cluster counts, gyration radii, percolation flags) are hand-computable. See the
oracle notes on each fixture.
"""

import os

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
from nexus.core.cluster import Cluster
from nexus.core.node import Node
from nexus.io.reader.xyz_reader import XYZReader

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


# --------------------------------------------------------------------------- #
# Paths to fixture trajectories
# --------------------------------------------------------------------------- #
@pytest.fixture
def data_dir():
    """Absolute path to ``tests/data/``."""
    return DATA_DIR


def _data_path(name):
    return os.path.join(DATA_DIR, name)


@pytest.fixture
def two_atoms_path():
    """2 type-A atoms, separation 1.0, box 10. Bonded at cutoff 1.1 -> 1 cluster."""
    return _data_path("two_atoms.xyz")


@pytest.fixture
def linear_chain_path():
    """5 colinear type-A atoms spaced 1.0, box 10. One cluster of size 5."""
    return _data_path("linear_chain.xyz")


@pytest.fixture
def two_clusters_path():
    """Two disjoint groups of 3 atoms, box 20. Exactly 2 clusters of size 3."""
    return _data_path("two_clusters.xyz")


@pytest.fixture
def pbc_pair_path():
    """2 atoms at x=0.25 and x=9.75, box 10. MIC separation 0.5 -> bonded only with PBC."""
    return _data_path("pbc_pair.xyz")


@pytest.fixture
def percolating_path():
    """27 atoms on a 3x3x3 simple-cubic lattice, box 3.0. With PBC at cutoff 1.1:
    one cluster percolating in x, y and z."""
    return _data_path("percolating_3x3x3.xyz")


@pytest.fixture
def two_frames_path():
    """Two-frame trajectory (2 atoms each). For scan/iteration tests."""
    return _data_path("two_frames.xyz")


@pytest.fixture
def malformed_path():
    """Single frame whose header lacks a Lattice="..." string. For reader error tests."""
    return _data_path("malformed.xyz")


@pytest.fixture
def sio2_chain_path():
    """Si-O bridged chain: Si-O-Si-O-Si-O-Si (4 Si, 3 bridging O), box 20.
    Si-O spacing 1.5. With Si-O cutoff ~2.0 the bond criterion (Si-O-Si) links all
    4 Si into one cluster; interior Si have O-coordination 2, terminal Si have 1."""
    return _data_path("sio2_chain.xyz")


@pytest.fixture
def lammps_path():
    """Two-frame LAMMPS dump (.lammpstrj), 2 type-A atoms each, box 10."""
    return _data_path("lammps_two.lammpstrj")


@pytest.fixture
def sio2_384_path():
    """5 frames of 384-atom amorphous SiO2 at 4.109 g/cc (high compression),
    from the thesis dataset (pos140B). With the bond Si-O-Si / coordination-O
    [4,6] alternating config, Si is mostly 5/6-coordinated and the SiO_5 / SiO_6
    networks percolate while sparse SiO_4 does not. Used for analyzer validation
    against golden values on real data."""
    return _data_path("sio2_384_dens4.109.xyz")


# --------------------------------------------------------------------------- #
# Settings factories
# --------------------------------------------------------------------------- #
def make_distance_settings(
    file_location,
    node_types=("A",),
    node_masses=None,
    cutoff=1.1,
    apply_pbc=True,
    with_all_analyzers=True,
    export_directory="test_export",
):
    """Build a validated ``Settings`` for single-/multi-type distance clustering.

    Mirrors how the quickstart scripts construct settings via ``SettingsBuilder``.
    Connectivity is ``[type, type]`` (distance criterion needs exactly 2 elements);
    one cutoff is added per type pair using the same ``cutoff`` value.
    """
    node_types = list(node_types)
    if node_masses is None:
        node_masses = [1.0] * len(node_types)

    cutoffs = []
    for i, t1 in enumerate(node_types):
        for t2 in node_types[i:]:
            cutoffs.append(Cutoff(type1=t1, type2=t2, distance=cutoff))

    connectivity = (
        [node_types[0], node_types[0]]
        if len(node_types) == 1
        else [node_types[0], node_types[-1]]
    )

    general = GeneralSettings(
        project_name="test",
        export_directory=export_directory,
        file_location=file_location,
        range_of_frames=(0, -1),
        apply_pbc=apply_pbc,
        verbose=False,
        save_logs=False,
        save_performance=False,
    )
    clustering = ClusteringSettings(
        criterion="distance",
        node_types=node_types,
        node_masses=node_masses,
        connectivity=connectivity,
        cutoffs=cutoffs,
    )
    analysis = AnalysisSettings(with_all=True)

    return (
        SettingsBuilder()
        .with_general(general)
        .with_lattice(LatticeSettings(apply_custom_lattice=False))
        .with_clustering(clustering)
        .with_analysis(analysis)
        .build()
    )


@pytest.fixture
def settings_factory():
    """Return the ``make_distance_settings`` factory for tests needing custom config."""
    return make_distance_settings


@pytest.fixture
def default_settings(two_atoms_path):
    """A ready distance-clustering ``Settings`` pointed at the two-atom fixture."""
    return make_distance_settings(two_atoms_path)


# --------------------------------------------------------------------------- #
# Frame factories
# --------------------------------------------------------------------------- #
def load_frame(file_location, settings, frame_id=0):
    """Parse one frame from an xyz fixture and initialize its nodes.

    Replicates the real pipeline: ``XYZReader.parse`` yields a raw frame, then
    ``initialize_nodes`` builds the filtered node list.
    """
    reader = XYZReader(settings)
    reader.scan()
    frame = next(reader.parse(frame_id))
    frame.initialize_nodes()
    return frame


@pytest.fixture
def frame_factory():
    """Return ``load_frame`` for tests that build frames from arbitrary fixtures."""
    return load_frame


@pytest.fixture
def two_atom_frame(two_atoms_path):
    """Frame with 2 bonded type-A atoms (PBC distance config)."""
    settings = make_distance_settings(two_atoms_path)
    return load_frame(two_atoms_path, settings)


@pytest.fixture
def linear_chain_frame(linear_chain_path):
    """Frame with 5 colinear type-A atoms."""
    settings = make_distance_settings(linear_chain_path)
    return load_frame(linear_chain_path, settings)


@pytest.fixture
def two_cluster_frame(two_clusters_path):
    """Frame with two disjoint 3-atom groups, box 20."""
    settings = make_distance_settings(two_clusters_path)
    return load_frame(two_clusters_path, settings)


@pytest.fixture
def percolating_frame(percolating_path):
    """Frame with a 3x3x3 simple-cubic lattice that percolates under PBC."""
    settings = make_distance_settings(percolating_path)
    return load_frame(percolating_path, settings)


# --------------------------------------------------------------------------- #
# Low-level object builders (no I/O)
# --------------------------------------------------------------------------- #
@pytest.fixture
def cubic_lattice():
    """A 10x10x10 orthorhombic lattice matrix."""
    return np.diag([10.0, 10.0, 10.0]).astype(float)


@pytest.fixture
def make_node():
    """Factory: ``make_node(symbol="A", node_id=0, position=(x,y,z), mass=1.0)``."""

    def _make(symbol="A", node_id=0, position=(0.0, 0.0, 0.0), mass=1.0):
        return Node(
            symbol=symbol,
            node_id=node_id,
            position=np.asarray(position, dtype=float),
            mass=mass,
        )

    return _make


@pytest.fixture
def make_cluster(cubic_lattice, two_atoms_path):
    """Factory: build a ``Cluster`` from a list of ``Node`` objects.

    ``make_cluster(nodes, connectivity="A-A", lattice=None)`` -> populated Cluster
    with ``size`` and ``root_id`` derived from the supplied nodes.
    """

    def _make(nodes, connectivity="A-A", lattice=None):
        lattice = cubic_lattice if lattice is None else lattice
        root_id = nodes[0].node_id if nodes else 0
        settings = make_distance_settings(two_atoms_path)
        cluster = Cluster(
            connectivity=connectivity,
            root_id=root_id,
            size=len(nodes),
            settings=settings,
            lattice=lattice,
        )
        for node in nodes:
            cluster.add_node(node)
        return cluster

    return _make
