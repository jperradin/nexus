"""Unit tests for nexus.core.cluster.Cluster.

Methods needing only cluster state (size, concentration, order parameter,
gyration radius from preset positions, MIC unwrap) are tested directly. The
BFS unwrapping is exercised on a hand-wired 2-node periodic graph so the
continuous-coordinate reconstruction has a known oracle.
"""

import numpy as np


class TestClusterBasics:
    def test_add_node_assigns_cluster_id(self, make_node, make_cluster):
        nodes = [make_node(node_id=0), make_node(node_id=1)]
        cluster = make_cluster(nodes)
        # root_id == first node's id; add_node tags every member with it.
        assert cluster.root_id == 0
        assert all(n.cluster_id == 0 for n in cluster.get_nodes())

    def test_get_size_matches_node_count(self, make_node, make_cluster):
        nodes = [make_node(node_id=i) for i in range(3)]
        assert make_cluster(nodes).get_size() == 3

    def test_get_connectivity(self, make_node, make_cluster):
        cluster = make_cluster([make_node(node_id=0)], connectivity="Si-Si")
        assert cluster.get_connectivity() == "Si-Si"


class TestConcentration:
    def test_ratio_size_over_total(self, make_node, make_cluster):
        cluster = make_cluster([make_node(node_id=i) for i in range(4)])
        cluster.total_nodes = 16
        cluster.calculate_concentration()
        assert cluster.concentration == 0.25

    def test_zero_total_leaves_default(self, make_node, make_cluster):
        cluster = make_cluster([make_node(node_id=0), make_node(node_id=1)])
        cluster.total_nodes = 0
        cluster.calculate_concentration()
        assert cluster.concentration == 0.0


class TestOrderParameter:
    def _cluster(self, make_node, make_cluster, n=4, total=8):
        cluster = make_cluster([make_node(node_id=i) for i in range(n)])
        cluster.total_nodes = total
        return cluster

    def test_one_direction(self, make_node, make_cluster):
        c = self._cluster(make_node, make_cluster)  # p_inf = 4/8 = 0.5
        c.percolation_probability = "x"
        c.calculate_order_parameter()
        assert c.order_parameter == [0.5, 0.0, 0.0]

    def test_two_directions(self, make_node, make_cluster):
        c = self._cluster(make_node, make_cluster)
        c.percolation_probability = "xy"
        c.calculate_order_parameter()
        assert c.order_parameter == [0.5, 0.5, 0.0]

    def test_three_directions(self, make_node, make_cluster):
        c = self._cluster(make_node, make_cluster)
        c.percolation_probability = "xyz"
        c.calculate_order_parameter()
        assert c.order_parameter == [0.5, 0.5, 0.5]

    def test_zero_total_leaves_default(self, make_node, make_cluster):
        c = self._cluster(make_node, make_cluster, total=0)
        c.percolation_probability = "xyz"
        c.calculate_order_parameter()
        assert c.order_parameter == [0.0, 0.0, 0.0]


class TestGyrationRadius:
    def test_single_node_is_zero(self, make_node, make_cluster):
        cluster = make_cluster([make_node(node_id=0)])
        cluster.calculate_gyration_radius()
        assert cluster.gyration_radius == 0.0

    def test_two_symmetric_nodes(self, make_node, make_cluster):
        # Unwrapped at (4,5,5) and (6,5,5), box 10 -> COM (5,5,5), no wrap shift.
        # Rg = rms distance from COM = sqrt((1^2 + 1^2)/2) = 1.0
        cluster = make_cluster([make_node(node_id=0), make_node(node_id=1)])
        cluster.unwrapped_positions = np.array([[4.0, 5.0, 5.0], [6.0, 5.0, 5.0]])
        cluster.calculate_gyration_radius()
        assert np.isclose(cluster.gyration_radius, 1.0)

    def test_empty_positions_is_zero(self, make_node, make_cluster):
        cluster = make_cluster([make_node(node_id=0), make_node(node_id=1)])
        cluster.unwrapped_positions = np.array([])
        cluster.calculate_gyration_radius()
        assert cluster.gyration_radius == 0.0


class TestCenterOfMass:
    def test_mean_inside_box_unshifted(self, make_node, make_cluster):
        cluster = make_cluster([make_node(node_id=0), make_node(node_id=1)])
        cluster.unwrapped_positions = np.array([[4.0, 5.0, 5.0], [6.0, 5.0, 5.0]])
        cluster.calculate_center_of_mass()
        np.testing.assert_allclose(cluster.center_of_mass, [5.0, 5.0, 5.0])

    def test_empty_positions_zero_com(self, make_node, make_cluster):
        cluster = make_cluster([make_node(node_id=0)])
        cluster.unwrapped_positions = np.array([])
        cluster.calculate_center_of_mass()
        np.testing.assert_allclose(cluster.center_of_mass, [0.0, 0.0, 0.0])


class TestUnwrapVector:
    """Minimum-image convention via fractional round."""

    def test_long_vector_folds_to_short(self, make_node, make_cluster):
        cluster = make_cluster([make_node(node_id=0)])  # box 10
        np.testing.assert_allclose(
            cluster._unwrap_vector(np.array([9.0, 0.0, 0.0])), [-1.0, 0.0, 0.0]
        )

    def test_short_vector_unchanged(self, make_node, make_cluster):
        cluster = make_cluster([make_node(node_id=0)])
        np.testing.assert_allclose(
            cluster._unwrap_vector(np.array([1.0, 0.0, 0.0])), [1.0, 0.0, 0.0]
        )

    def test_negative_long_vector(self, make_node, make_cluster):
        cluster = make_cluster([make_node(node_id=0)])
        np.testing.assert_allclose(
            cluster._unwrap_vector(np.array([-9.0, 0.0, 0.0])), [1.0, 0.0, 0.0]
        )


class TestUnwrappedPositionsBFS:
    def test_two_node_periodic_pair(self, make_node, make_cluster):
        # Two atoms bonded across the x boundary: x=0.25 and x=9.75, box 10.
        # MIC separation = 0.5, so unwrapping must place them 0.5 apart with no
        # net period vector (a simple 2-node chain does not loop).
        n0 = make_node(node_id=0, position=(0.25, 5.0, 5.0))
        n1 = make_node(node_id=1, position=(9.75, 5.0, 5.0))
        n0.add_neighbor(n1)
        n1.add_neighbor(n0)
        cluster = make_cluster([n0, n1])  # criterion = "distance"

        cluster.calculate_unwrapped_positions()

        assert cluster.linkages == [(0, 1)]
        assert len(cluster.period_vectors) == 0
        sep = np.linalg.norm(
            cluster.unwrapped_positions[0] - cluster.unwrapped_positions[1]
        )
        assert np.isclose(sep, 0.5)


class TestPercolationMath:
    def test_single_node_not_percolating(self, make_node, make_cluster):
        cluster = make_cluster([make_node(node_id=0)])
        cluster.size = 1
        cluster.calculate_percolation_probability()
        assert cluster.is_percolating is False
        assert cluster.percolation_probability == ""

    def test_period_dimension_rank_one(self, make_node, make_cluster):
        cluster = make_cluster([make_node(node_id=0)])
        cluster.period_vectors = [np.array([10.0, 0.0, 0.0])]
        assert cluster._calculate_period_dimension() == 1

    def test_period_dimension_rank_three(self, make_node, make_cluster):
        cluster = make_cluster([make_node(node_id=0)])
        cluster.period_vectors = [
            np.array([10.0, 0.0, 0.0]),
            np.array([0.0, 10.0, 0.0]),
            np.array([0.0, 0.0, 10.0]),
        ]
        assert cluster._calculate_period_dimension() == 3

    def test_directions_all_three(self, make_node, make_cluster):
        cluster = make_cluster([make_node(node_id=0)])  # box 10
        cluster.period_vectors = [
            np.array([10.0, 0.0, 0.0]),
            np.array([0.0, 10.0, 0.0]),
            np.array([0.0, 0.0, 10.0]),
        ]
        assert cluster._get_percolation_directions() == "xyz"

    def test_directions_single_axis(self, make_node, make_cluster):
        cluster = make_cluster([make_node(node_id=0)])
        cluster.period_vectors = [np.array([10.0, 0.0, 0.0])]
        assert cluster._get_percolation_directions() == "x"

    def test_no_period_vectors(self, make_node, make_cluster):
        cluster = make_cluster([make_node(node_id=0)])
        cluster.period_vectors = []
        assert cluster._calculate_period_dimension() == 0
        assert cluster._get_percolation_directions() == ""
