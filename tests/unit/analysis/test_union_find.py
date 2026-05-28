"""Unit tests for the union-find core in BaseClusteringStrategy (find / union).

BaseClusteringStrategy is abstract, so a trivial concrete subclass is defined
here. find()/union() operate purely on Node.parent links and do not touch the
frame, so any initialized frame fixture suffices for construction.

Nodes must have distinct node_ids: Node is an order=True dataclass comparing on
(symbol, node_id, mass, coordination), and find() relies on ``parent != node``.
"""

from nexus.analysis.strategies.base_strategy import BaseClusteringStrategy


class _UF(BaseClusteringStrategy):
    """Minimal concrete strategy exposing find()/union()."""

    def build_clusters(self):
        return []


def make_uf(two_atom_frame, default_settings):
    return _UF(two_atom_frame, default_settings)


class TestFind:
    def test_isolated_node_is_own_root(self, two_atom_frame, default_settings, make_node):
        uf = make_uf(two_atom_frame, default_settings)
        n = make_node(node_id=0)
        assert uf.find(n) is n

    def test_path_compression(self, two_atom_frame, default_settings, make_node):
        uf = make_uf(two_atom_frame, default_settings)
        a = make_node(node_id=0)
        b = make_node(node_id=1)
        c = make_node(node_id=2)
        a.parent = b
        b.parent = c  # chain a -> b -> c
        root = uf.find(a)
        assert root is c
        # After find, a points directly at the root (compressed).
        assert a.parent is c


class TestUnion:
    def test_union_merges_two(self, two_atom_frame, default_settings, make_node):
        uf = make_uf(two_atom_frame, default_settings)
        a = make_node(node_id=0)
        b = make_node(node_id=1)
        uf.union(a, b)
        assert uf.find(a) is uf.find(b)

    def test_union_first_arg_becomes_root(self, two_atom_frame, default_settings, make_node):
        uf = make_uf(two_atom_frame, default_settings)
        a = make_node(node_id=0)
        b = make_node(node_id=1)
        uf.union(a, b)  # root_2.parent = root_1 -> a is root
        assert uf.find(b) is a

    def test_transitive_merge(self, two_atom_frame, default_settings, make_node):
        uf = make_uf(two_atom_frame, default_settings)
        a = make_node(node_id=0)
        b = make_node(node_id=1)
        c = make_node(node_id=2)
        uf.union(a, b)
        uf.union(b, c)
        root = uf.find(a)
        assert uf.find(b) is root and uf.find(c) is root

    def test_union_idempotent(self, two_atom_frame, default_settings, make_node):
        uf = make_uf(two_atom_frame, default_settings)
        a = make_node(node_id=0)
        b = make_node(node_id=1)
        uf.union(a, b)
        uf.union(a, b)  # repeat must not create a cycle
        assert uf.find(a) is a
        assert uf.find(b) is a

    def test_union_same_set_no_change(self, two_atom_frame, default_settings, make_node):
        uf = make_uf(two_atom_frame, default_settings)
        a = make_node(node_id=0)
        b = make_node(node_id=1)
        uf.union(a, b)
        root_before = uf.find(a)
        uf.union(b, a)  # already together
        assert uf.find(a) is root_before
