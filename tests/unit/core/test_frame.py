"""Unit tests for nexus.core.frame.Frame."""

import numpy as np
import pytest

from nexus.core.frame import Frame


def raw_frame(symbols, positions, settings, lattice=None):
    """Build an un-initialized Frame straight from raw symbol/position data."""
    lattice = np.diag([10.0, 10.0, 10.0]).astype(float) if lattice is None else lattice
    data = {"symbol": symbols, "position": [np.asarray(p, dtype=float) for p in positions]}
    return Frame(
        frame_id=0,
        nodes=[],
        lattice=lattice,
        _data=data,
        _settings=settings,
    )


class TestInitializeNodes:
    def test_keeps_only_selected_species(self, default_settings):
        # default_settings.node_types == ["A"]; B atoms must be filtered out.
        frame = raw_frame(
            ["A", "B", "A"],
            [(0, 0, 0), (1, 0, 0), (2, 0, 0)],
            default_settings,
        )
        frame.initialize_nodes()
        assert len(frame) == 2
        assert all(n.symbol == "A" for n in frame.get_nodes())

    def test_assigns_sequential_ids(self, default_settings):
        frame = raw_frame(["A", "A", "A"], [(0, 0, 0), (1, 0, 0), (2, 0, 0)], default_settings)
        frame.initialize_nodes()
        assert [n.node_id for n in frame.get_nodes()] == [0, 1, 2]

    def test_mismatched_lengths_raise(self, default_settings):
        frame = raw_frame(["A", "A"], [(0, 0, 0)], default_settings)
        with pytest.raises(ValueError):
            frame.initialize_nodes()


class TestLattice:
    def test_set_and_get(self, two_atom_frame):
        new_lat = np.diag([5.0, 6.0, 7.0]).astype(float)
        two_atom_frame.set_lattice(new_lat)
        np.testing.assert_allclose(two_atom_frame.get_lattice(), new_lat)

    def test_wrong_shape_raises(self, two_atom_frame):
        with pytest.raises(ValueError):
            two_atom_frame.set_lattice(np.eye(2))

    def test_singular_lattice_raises(self, two_atom_frame):
        with pytest.raises(ValueError):
            two_atom_frame.set_lattice(np.zeros((3, 3)))


class TestNodeAccess:
    def test_unique_elements(self, two_atom_frame):
        assert list(two_atom_frame.get_unique_elements()) == ["A"]

    def test_get_node_by_id_found(self, linear_chain_frame):
        node = linear_chain_frame.get_node_by_id(3)
        assert node is not None and node.node_id == 3

    def test_get_node_by_id_missing_returns_none(self, two_atom_frame):
        assert two_atom_frame.get_node_by_id(999) is None

    def test_positions_shape(self, linear_chain_frame):
        pos = linear_chain_frame.get_positions()
        assert pos.shape == (5, 3)

    def test_wrapped_positions_inside_box(self, percolating_frame):
        wrapped = percolating_frame.get_wrapped_positions()
        diag = np.diag(percolating_frame.get_lattice())
        assert np.all(wrapped >= 0.0)
        assert np.all(wrapped < diag)

    def test_len(self, linear_chain_frame):
        assert len(linear_chain_frame) == 5


class TestClustersAndConnectivity:
    def test_add_cluster_initializes_list(self, two_atom_frame, make_cluster):
        cluster = make_cluster(two_atom_frame.get_nodes())
        two_atom_frame.clusters = None
        two_atom_frame.add_cluster(cluster)
        assert two_atom_frame.get_clusters() == [cluster]

    def test_networking_nodes_sums_cluster_sizes(self, two_atom_frame, make_cluster):
        cluster = make_cluster(two_atom_frame.get_nodes(), connectivity="A-A")
        two_atom_frame.set_clusters([cluster])
        assert two_atom_frame.get_networking_nodes() == 2

    def test_set_and_get_connectivities(self, two_atom_frame):
        two_atom_frame.set_connectivities(["A-A"])
        assert two_atom_frame.get_connectivities() == ["A-A"]

    def test_concentration_ratio(self, two_atom_frame, make_cluster):
        cluster = make_cluster(two_atom_frame.get_nodes(), connectivity="A-A")
        cluster.total_nodes = 2
        two_atom_frame.set_clusters([cluster])
        two_atom_frame.set_connectivities(["A-A"])
        conc = two_atom_frame.get_concentration()
        # 2 networking nodes (both type A) all in the cluster -> ratio 1.0
        assert conc["A-A"] == 1.0

    def test_concentration_zero_for_missing_connectivity(self, two_atom_frame, make_cluster):
        cluster = make_cluster(two_atom_frame.get_nodes(), connectivity="A-A")
        cluster.total_nodes = 2
        two_atom_frame.set_clusters([cluster])
        two_atom_frame.set_connectivities(["A-A", "B-B"])
        conc = two_atom_frame.get_concentration()
        assert conc["B-B"] == 0.0
