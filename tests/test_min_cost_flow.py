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
import gurobi_optimods.min_cost_flow as mcf


class TestMinCostFlow(unittest.TestCase):
    def test_min_cost_flow(self):
        edge_data, node_data = datasets.load_min_cost_flow()
        cost, sol = mcf.min_cost_flow(edge_data, node_data)
        self.assertEqual(cost, 31)
        self.assertEqual(sol.tolist(), [1.0, 1.0, 1.0, 0.0, 2.0, 0.0, 2.0])
        self.assertEqual(
            sol.index.tolist(), [(0, 1), (0, 2), (1, 3), (2, 3), (2, 4), (3, 5), (4, 5)]
        )

    def test_min_cost_flow_scipy(self):
        G, cap, cost, demands = datasets.load_min_cost_flow_scipy()
        cost, sol = mcf.min_cost_flow_scipy(G, cap, cost, demands)
        self.assertEqual(cost, 31)
        expected = np.array(
            [
                [0.0, 1.0, 1.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 2.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 2.0],
            ]
        )
        assert_array_equal(sol.toarray(), expected)

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_min_cost_flow_networkx(self):
        G = datasets.load_min_cost_flow_networkx()
        cost, sol = mcf.min_cost_flow_networkx(G)
        self.assertEqual(cost, 31)
        self.assertEqual(
            sol, {(0, 1): 1.0, (0, 2): 1.0, (1, 3): 1.0, (2, 4): 2.0, (4, 5): 2.0}
        )


class TestMinCostFlow2(unittest.TestCase):
    def test_min_cost_flow(self):
        edge_data, node_data = datasets.load_min_cost_flow2()
        cost, sol = mcf.min_cost_flow(edge_data, node_data)
        self.assertEqual(cost, 150)
        self.assertEqual(sol.tolist(), [12.0, 8.0, 4.0, 8.0, 0.0, 11.0, 5.0, 10.0, 0.0])
        self.assertEqual(
            sol.index.tolist(),
            [(0, 1), (0, 2), (1, 3), (1, 2), (1, 4), (2, 3), (2, 4), (3, 4), (4, 2)],
        )

    def test_min_cost_flow_scipy(self):
        G, cap, cost, demands = datasets.load_min_cost_flow2_scipy()
        cost, sol = mcf.min_cost_flow_scipy(G, cap, cost, demands)
        self.assertEqual(cost, 150)
        expected = np.array(
            [
                [0.0, 12.0, 8.0, 0.0, 0.0],
                [0.0, 0.0, 8.0, 4.0, 0.0],
                [0.0, 0.0, 0.0, 11.0, 5.0],
                [0.0, 0.0, 0.0, 0.0, 10.0],
                [0.0, 0.0, 0.0, 0.0, 0.0],
            ]
        )
        assert_array_equal(sol.toarray(), expected)

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_min_cost_flow_networkx(self):
        G = datasets.load_min_cost_flow2_networkx()
        cost, sol = mcf.min_cost_flow_networkx(G)
        self.assertEqual(cost, 150)
        self.assertEqual(
            sol,
            {
                (0, 1): 12.0,
                (0, 2): 8.0,
                (1, 3): 4.0,
                (1, 2): 8.0,
                (2, 3): 11.0,
                (2, 4): 5.0,
                (3, 4): 10.0,
            },
        )
