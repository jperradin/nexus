"""Unit tests for nexus.core.node.Node."""

import numpy as np

from nexus.core.node import Node


class TestNodeDefaults:
    """__post_init__ default-filling behavior."""

    def test_self_parent_by_default(self, make_node):
        # Each node starts as the root of its own union-find set.
        n = make_node(node_id=0)
        assert n.parent is n

    def test_cluster_id_defaults_to_node_id(self, make_node):
        n = make_node(node_id=7)
        assert n.cluster_id == 7

    def test_mass_defaults_to_one_when_none(self):
        n = Node(symbol="A", node_id=0, position=np.zeros(3), mass=None)
        assert n.mass == 1.0

    def test_coordination_defaults_to_zero(self):
        n = Node(symbol="A", node_id=0, position=np.zeros(3))
        assert n.coordination == 0

    def test_neighbors_and_other_default_empty(self):
        n = Node(symbol="A", node_id=0, position=np.zeros(3))
        assert n.neighbors == []
        assert n.other == []


class TestAddNeighbor:
    def test_appends_neighbor(self, make_node):
        a = make_node(node_id=0)
        b = make_node(node_id=1, position=(1.0, 0.0, 0.0))
        a.add_neighbor(b)
        assert a.neighbors == [b]

    def test_multiple_neighbors_preserve_order(self, make_node):
        a = make_node(node_id=0)
        b = make_node(node_id=1)
        c = make_node(node_id=2)
        a.add_neighbor(b)
        a.add_neighbor(c)
        assert a.neighbors == [b, c]


class TestResetParent:
    def test_reset_restores_self_parent(self, make_node):
        a = make_node(node_id=0)
        b = make_node(node_id=1)
        a.parent = b  # simulate a union
        assert a.parent is b
        a.reset_parent()
        assert a.parent is a


class TestSetCoordination:
    def test_sets_value(self, make_node):
        n = make_node(node_id=0)
        n.set_coordination(4)
        assert n.coordination == 4


class TestWrapPosition:
    """Static wrap_position uses fractional coords + floor -> result in [0, L)."""

    def test_point_inside_box_unchanged(self, cubic_lattice):
        pos = np.array([1.0, 2.0, 3.0])
        np.testing.assert_allclose(Node.wrap_position(pos, cubic_lattice), pos)

    def test_point_above_box_wraps_in(self, cubic_lattice):
        # box side 10: x=12 -> 2, y=10 -> 0
        wrapped = Node.wrap_position(np.array([12.0, 10.0, 5.0]), cubic_lattice)
        assert np.all(wrapped >= 0.0) and np.all(wrapped < 10.0)
        np.testing.assert_allclose(wrapped, [2.0, 0.0, 5.0])

    def test_negative_coord_wraps_positive(self, cubic_lattice):
        wrapped = Node.wrap_position(np.array([-1.0, -3.0, 0.0]), cubic_lattice)
        np.testing.assert_allclose(wrapped, [9.0, 7.0, 0.0])


class TestReprAndOrdering:
    def test_str_contains_id_and_symbol(self, make_node):
        s = str(make_node(symbol="Si", node_id=3))
        assert "Si" in s and "3" in s

    def test_ordering_by_compare_fields(self, make_node):
        # @dataclass(order=True): compare tuple = (symbol, node_id, mass, coordination).
        a = make_node(symbol="A", node_id=0)
        b = make_node(symbol="A", node_id=1)
        assert a < b
