import unittest

import pandas as pd

try:
    import networkx as nx
except ImportError:
    nx = None

import gurobi_optimods.datasets as datasets
from gurobi_optimods.min_cut import min_cut
from gurobi_optimods.datasets import (
    _convert_pandas_to_digraph,
    _convert_pandas_to_scipy,
)


class TestMinCut(unittest.TestCase):
    def setUp(self):
        self.expected_partition = (set({0, 1}), set({2, 3, 4, 5}))
        self.expected_cut_value = 3.0
        self.expected_cut_set = set([(0, 2), (1, 3)])

    def test_infeasible(self):
        edge_data, _ = datasets.load_graph()
        edge_data["capacity"] = [-1] * len(edge_data)
        with self.assertRaisesRegex(ValueError, "Unsatisfiable flows"):
            obj, sol, cutset = min_cut(edge_data, 0, 5)

    def test_pandas(self):
        edge_data, node_data = datasets.load_graph()
        cut_value, partition, cutset = min_cut(edge_data, 0, 5)
        self.assertEqual(cut_value, self.expected_cut_value)
        self.assertEqual(partition[0], self.expected_partition[0])
        self.assertEqual(partition[1], self.expected_partition[1])
        self.assertEqual(cutset, self.expected_cut_set)

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_networkx(self):
        G = datasets.load_graph_networkx()
        cut_value, partition, cutset = min_cut(G, 0, 5)
        self.assertEqual(cut_value, self.expected_cut_value)
        self.assertEqual(partition[0], self.expected_partition[0])
        self.assertEqual(partition[1], self.expected_partition[1])
        self.assertEqual(cutset, self.expected_cut_set)

    def test_scipy(self):
        G, capacity, _, _ = datasets.load_graph_scipy()
        G.data = capacity.data
        cut_value, partition, cutset = min_cut(G, 0, 5)
        self.assertEqual(cut_value, self.expected_cut_value)
        self.assertEqual(partition[0], self.expected_partition[0])
        self.assertEqual(partition[1], self.expected_partition[1])
        self.assertEqual(cutset, self.expected_cut_set)


class TestMinCut2(unittest.TestCase):
    def setUp(self):
        self.expected_partition = (set({0}), set({1, 2, 3, 4}))
        self.expected_cut_value = 23.0
        self.expected_cut_set = set([(0, 1), (0, 2)])

    def test_pandas(self):
        edge_data, _ = datasets.load_graph2()
        cut_value, partition, cutset = min_cut(edge_data, 0, 4)
        self.assertEqual(cut_value, self.expected_cut_value)
        self.assertEqual(partition[0], self.expected_partition[0])
        self.assertEqual(partition[1], self.expected_partition[1])
        self.assertEqual(cutset, self.expected_cut_set)

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_networkx(self):
        G = datasets.load_graph2_networkx()
        cut_value, partition, cutset = min_cut(G, 0, 4)
        self.assertEqual(cut_value, self.expected_cut_value)
        self.assertEqual(partition[0], self.expected_partition[0])
        self.assertEqual(partition[1], self.expected_partition[1])
        self.assertEqual(cutset, self.expected_cut_set)

    def test_scipy(self):
        G, capacity, _, _ = datasets.load_graph2_scipy()
        G.data = capacity.data
        cut_value, partition, cutset = min_cut(G, 0, 4)
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
        cut_value, partition, cutset = min_cut(self.arc_data, 0, 6)
        self.assertEqual(cut_value, self.expected_cut_value)
        self.assertEqual(partition[0], self.expected_partition[0])
        self.assertEqual(partition[1], self.expected_partition[1])
        self.assertEqual(cutset, self.expected_cut_set)

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_networkx(self):
        cut_value, partition, cutset = min_cut(
            _convert_pandas_to_digraph(self.arc_data, None, demand=False), 0, 6
        )
        self.assertEqual(cut_value, self.expected_cut_value)
        self.assertEqual(partition[0], self.expected_partition[0])
        self.assertEqual(partition[1], self.expected_partition[1])
        self.assertEqual(cutset, self.expected_cut_set)

    def test_scipy(self):
        G, capacity, _, _ = _convert_pandas_to_scipy(
            self.arc_data, None, cost=False, demand=False
        )

        G.data = capacity.data
        cut_value, partition, cutset = min_cut(G, 0, 6)
        self.assertEqual(cut_value, self.expected_cut_value)
        self.assertEqual(partition[0], self.expected_partition[0])
        self.assertEqual(partition[1], self.expected_partition[1])
        self.assertEqual(cutset, self.expected_cut_set)
