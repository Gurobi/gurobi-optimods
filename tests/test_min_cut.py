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
