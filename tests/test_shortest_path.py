import unittest

import numpy as np
import pandas as pd
import scipy.sparse as sp

try:
    import networkx as nx
except ImportError:
    nx = None

from numpy.testing import assert_array_equal

import gurobi_optimods.datasets as datasets
import gurobi_optimods.shortest_path as mcf


class TestShortestPath(unittest.TestCase):
    def test_pandas(self):
        edge_data, _ = datasets.load_min_cost_flow()
        cost, sol = mcf.shortest_path(edge_data, 0, 5)
        self.assertEqual(cost, 11)
        sol = sol[sol > 0]
        self.assertEqual(sol.tolist(), [1.0, 1.0, 1.0])
        self.assertEqual(sol.index.tolist(), [(0, 1), (1, 3), (3, 5)])

    def test_infeasible(self):
        G, _, _, _ = datasets.load_min_cost_flow_scipy()
        with self.assertRaisesRegex(ValueError, "Unsatisfiable flows"):
            obj, sol = mcf.shortest_path(G, 0, 0)

    def test_scipy(self):
        G, _, cost, _ = datasets.load_min_cost_flow_scipy()
        G.data = cost.data
        cost, sol = mcf.shortest_path(G, 0, 5)
        self.assertEqual(cost, 11.0)
        expected = np.array(
            [
                [0, 1, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 1],
                [0, 0, 0, 0, 0, 0],
            ]
        )
        assert_array_equal(sol.toarray(), expected)

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_networkx(self):
        G = datasets.load_min_cost_flow_networkx()
        m = nx.to_scipy_sparse_array(G).tocsr()
        cost, sol_graph = mcf.shortest_path(G, 0, 5)
        self.assertEqual(cost, 11)
        self.assertEqual(list(sol_graph.edges()), [(0, 1), (1, 3), (3, 5)])


class TestShortestPath2(unittest.TestCase):
    def test_pandas(self):
        edge_data, _ = datasets.load_min_cost_flow2()
        cost, sol = mcf.shortest_path(edge_data, 0, 4)
        sol = sol[sol > 0]
        self.assertEqual(cost, 7.0)
        self.assertEqual(sol.tolist(), [1.0] * 3)
        self.assertEqual(
            sol.index.tolist(),
            [(0, 2), (2, 3), (3, 4)],
        )

    def test_sp(self):
        G, _, cost, _ = datasets.load_min_cost_flow2_scipy()
        G.data = cost.data
        cost, sol = mcf.shortest_path(G, 0, 4)
        self.assertEqual(cost, 7.0)
        expected = np.array(
            [
                [0, 0, 1, 0, 0],
                [0, 0, 0, 0, 0],
                [0, 0, 0, 1, 0],
                [0, 0, 0, 0, 1],
                [0, 0, 0, 0, 0],
            ]
        )
        assert_array_equal(sol.toarray(), expected)

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_networkx(self):
        G = datasets.load_min_cost_flow2_networkx()
        cost, sol = mcf.shortest_path(G, 0, 4)
        self.assertEqual(cost, 7.0)
        self.assertEqual(list(sol.edges()), [(0, 2), (2, 3), (3, 4)])
