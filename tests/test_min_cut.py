import unittest

import numpy as np
import pandas as pd

try:
    import networkx as nx
except ImportError:
    nx = None

import gurobi_optimods.datasets as datasets
from gurobi_optimods.datasets import (
    _convert_pandas_to_digraph,
    _convert_pandas_to_scipy,
)
from gurobi_optimods.min_cut import min_cut

from .test_min_cost_flow import (
    load_graph2_networkx,
    load_graph2_pandas,
    load_graph2_scipy,
)


class TestMinCut(unittest.TestCase):
    def setUp(self):
        self.expected_partition = (set({0, 1}), set({2, 3, 4, 5}))
        self.expected_cut_value = 3.0
        self.expected_cut_set = set([(0, 2), (1, 3)])

    def test_pandas(self):
        edge_data, node_data = datasets.simple_graph_pandas()
        res = min_cut(edge_data, 0, 5)
        cut_value, partition, cutset = res.cut_value, res.partition, res.cutset
        self.assertEqual(cut_value, self.expected_cut_value)
        self.assertEqual(partition[0], self.expected_partition[0])
        self.assertEqual(partition[1], self.expected_partition[1])
        self.assertEqual(cutset, self.expected_cut_set)

    def test_empty_pandas(self):
        edge_data, _ = datasets.simple_graph_pandas()
        edge_data["capacity"] = [0] * len(edge_data)
        res = min_cut(edge_data, 0, 5)
        cut_value, partition, cutset = res.cut_value, res.partition, res.cutset
        self.assertEqual(cut_value, 0.0)
        self.assertEqual(partition, (set(), set()))
        self.assertEqual(cutset, set())

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_networkx(self):
        G = datasets.simple_graph_networkx()
        res = min_cut(G, 0, 5)
        cut_value, partition, cutset = res.cut_value, res.partition, res.cutset
        self.assertEqual(cut_value, self.expected_cut_value)
        self.assertEqual(partition[0], self.expected_partition[0])
        self.assertEqual(partition[1], self.expected_partition[1])
        self.assertEqual(cutset, self.expected_cut_set)

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_empty_networkx(self):
        G = datasets.simple_graph_networkx()
        nx.set_edge_attributes(G, 0, "capacity")
        res = min_cut(G, 0, 5)
        cut_value, partition, cutset = res.cut_value, res.partition, res.cutset
        self.assertEqual(cut_value, 0.0)
        self.assertEqual(partition, (set(), set()))
        self.assertEqual(cutset, set())

    def test_scipy(self):
        G, capacity, _, _ = datasets.simple_graph_scipy()
        G.data = capacity.data
        res = min_cut(G, 0, 5)
        cut_value, partition, cutset = res.cut_value, res.partition, res.cutset
        self.assertEqual(cut_value, self.expected_cut_value)
        self.assertEqual(partition[0], self.expected_partition[0])
        self.assertEqual(partition[1], self.expected_partition[1])
        self.assertEqual(cutset, self.expected_cut_set)

    def test_empty_scipy(self):
        G, capacity, _, _ = datasets.simple_graph_scipy()
        G.data = np.repeat(1e-20, len(G.data))
        res = min_cut(G, 0, 5)
        cut_value, partition, cutset = res.cut_value, res.partition, res.cutset
        self.assertEqual(cut_value, 0.0)
        self.assertEqual(partition, (set(), set()))
        self.assertEqual(cutset, set())


class TestMinCut2(unittest.TestCase):
    def setUp(self):
        self.expected_partition = (set({0}), set({1, 2, 3, 4}))
        self.expected_cut_value = 23.0
        self.expected_cut_set = set([(0, 1), (0, 2)])

    def test_pandas(self):
        edge_data, _ = load_graph2_pandas()
        res = min_cut(edge_data, 0, 4)
        cut_value, partition, cutset = res.cut_value, res.partition, res.cutset
        self.assertEqual(cut_value, self.expected_cut_value)
        self.assertEqual(partition[0], self.expected_partition[0])
        self.assertEqual(partition[1], self.expected_partition[1])
        self.assertEqual(cutset, self.expected_cut_set)

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_networkx(self):
        G = load_graph2_networkx()
        res = min_cut(G, 0, 4)
        cut_value, partition, cutset = res.cut_value, res.partition, res.cutset
        self.assertEqual(cut_value, self.expected_cut_value)
        self.assertEqual(partition[0], self.expected_partition[0])
        self.assertEqual(partition[1], self.expected_partition[1])
        self.assertEqual(cutset, self.expected_cut_set)

    def test_scipy(self):
        G, capacity, _, _ = load_graph2_scipy()
        G.data = capacity.data
        res = min_cut(G, 0, 4)
        cut_value, partition, cutset = res.cut_value, res.partition, res.cutset
        self.assertEqual(cut_value, self.expected_cut_value)
        self.assertEqual(partition[0], self.expected_partition[0])
        self.assertEqual(partition[1], self.expected_partition[1])
        self.assertEqual(cutset, self.expected_cut_set)


class TestMinCut3(unittest.TestCase):
    def setUp(self):
        self.expected_partition = (set({3, 1, 0}), set({4, 6, 5, 2}))
        self.expected_cut_value = 3.0
        self.expected_cut_set = set([(0, 2), (3, 6)])
        self.arc_data = pd.DataFrame(
            [
                {"source": 0, "target": 1, "capacity": 5.0},
                {"source": 0, "target": 2, "capacity": 1.0},
                {"source": 1, "target": 3, "capacity": 3.0},
                {"source": 2, "target": 3, "capacity": 5.0},
                {"source": 2, "target": 4, "capacity": 4.0},
                {"source": 4, "target": 5, "capacity": 2.0},
                {"source": 3, "target": 6, "capacity": 2.0},
                {"source": 5, "target": 6, "capacity": 3.0},
            ]
        ).set_index(["source", "target"])

    def test_pandas(self):
        res = min_cut(self.arc_data, 0, 6)
        cut_value, partition, cutset = res.cut_value, res.partition, res.cutset
        self.assertEqual(cut_value, self.expected_cut_value)
        self.assertEqual(partition[0], self.expected_partition[0])
        self.assertEqual(partition[1], self.expected_partition[1])
        self.assertEqual(cutset, self.expected_cut_set)

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_networkx(self):
        res = min_cut(
            _convert_pandas_to_digraph(self.arc_data, None, demand=False), 0, 6
        )
        cut_value, partition, cutset = res.cut_value, res.partition, res.cutset
        self.assertEqual(cut_value, self.expected_cut_value)
        self.assertEqual(partition[0], self.expected_partition[0])
        self.assertEqual(partition[1], self.expected_partition[1])
        self.assertEqual(cutset, self.expected_cut_set)

    def test_scipy(self):
        G, capacity, _, _ = _convert_pandas_to_scipy(
            self.arc_data, None, cost=False, demand=False
        )

        G.data = capacity.data
        res = min_cut(G, 0, 6)
        cut_value, partition, cutset = res.cut_value, res.partition, res.cutset
        self.assertEqual(cut_value, self.expected_cut_value)
        self.assertEqual(partition[0], self.expected_partition[0])
        self.assertEqual(partition[1], self.expected_partition[1])
        self.assertEqual(cutset, self.expected_cut_set)


class TestMinCutBindingDummyEdge(unittest.TestCase):
    """Single edge 0 -> 1 with capacity 10: the dummy edge capacity (the
    max flow bound) is binding and must not be reported as the cut."""

    def check(self, res):
        self.assertEqual(res.cut_value, 10.0)
        self.assertEqual(res.partition, ({0}, {1}))
        self.assertEqual(res.cutset, {(0, 1)})

    def test_scipy(self):
        from scipy import sparse

        G = sparse.lil_array((2, 2))
        G[0, 1] = 10
        self.check(min_cut(G.tocoo(), 0, 1))

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_networkx(self):
        G = nx.DiGraph()
        G.add_edge(0, 1, capacity=10)
        self.check(min_cut(G, 0, 1))

    def test_pandas(self):
        arc_data = pd.DataFrame(
            [{"source": 0, "target": 1, "capacity": 10.0}]
        ).set_index(["source", "target"])
        self.check(min_cut(arc_data, 0, 1))


class TestMinCutBottleneckAtTerminal(unittest.TestCase):
    """Non-trivial graphs where the min cut is exactly the set of edges at the
    source or sink, so the dummy edge capacity is binding."""

    # Min cut is the single source edge (0, 1), value 7. Source is node 0,
    # sink is node 4.
    source_cut_edges = [
        (0, 1, 7.0),
        (1, 2, 5.0),
        (1, 3, 2.0),
        (1, 4, 5.0),
        (2, 3, 2.0),
        (2, 4, 3.0),
        (3, 4, 5.0),
    ]

    # Min cut is the single sink edge (4, 5), value 5. Source is node 0, sink
    # is node 5. The source capacity (13) is larger than the sink capacity (5),
    # so only the pandas bound (min of the two) is binding.
    sink_cut_edges = [
        (0, 1, 6.0),
        (0, 2, 5.0),
        (0, 3, 1.0),
        (0, 4, 1.0),
        (1, 3, 9.0),
        (1, 4, 9.0),
        (2, 3, 2.0),
        (3, 4, 7.0),
        (4, 5, 5.0),
    ]

    def check_source_cut(self, res):
        self.assertEqual(res.cut_value, 7.0)
        self.assertEqual(res.partition, ({0}, {1, 2, 3, 4}))
        self.assertEqual(res.cutset, {(0, 1)})

    def check_sink_cut(self, res):
        self.assertEqual(res.cut_value, 5.0)
        self.assertEqual(res.partition, ({0, 1, 2, 3, 4}, {5}))
        self.assertEqual(res.cutset, {(4, 5)})

    def test_scipy_source_cut(self):
        from scipy import sparse

        rows, cols, capacities = zip(*self.source_cut_edges)
        G = sparse.coo_array((capacities, (rows, cols)), shape=(5, 5))
        self.check_source_cut(min_cut(G, 0, 4))

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_networkx_source_cut(self):
        G = nx.DiGraph()
        for i, j, capacity in self.source_cut_edges:
            G.add_edge(i, j, capacity=capacity)
        self.check_source_cut(min_cut(G, 0, 4))

    def test_pandas_source_cut(self):
        arc_data = pd.DataFrame(
            self.source_cut_edges, columns=["source", "target", "capacity"]
        ).set_index(["source", "target"])
        self.check_source_cut(min_cut(arc_data, 0, 4))

    def test_pandas_sink_cut(self):
        arc_data = pd.DataFrame(
            self.sink_cut_edges, columns=["source", "target", "capacity"]
        ).set_index(["source", "target"])
        self.check_sink_cut(min_cut(arc_data, 0, 5))


class TestMinCutScipyNonzeroSource(unittest.TestCase):
    """Source is node 1, and node 0 has no outgoing edges. The max flow bound
    must be computed from the source row, not row 0. The min cut (4) is below
    the bound (10), so the dummy edge is not binding."""

    def test_scipy(self):
        from scipy import sparse

        G = sparse.lil_array((4, 4))
        G[1, 2] = 10
        G[2, 3] = 4
        res = min_cut(G.tocoo(), 1, 3)
        self.assertEqual(res.cut_value, 4.0)
        self.assertEqual(res.partition, ({1, 2}, {0, 3}))
        self.assertEqual(res.cutset, {(2, 3)})
