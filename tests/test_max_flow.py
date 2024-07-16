import unittest

import numpy as np

try:
    import networkx as nx
except ImportError:
    nx = None

import pandas as pd

import gurobi_optimods.datasets as datasets
from gurobi_optimods.max_flow import max_flow

from .test_graph_utils import (
    check_solution_networkx,
    check_solution_pandas,
    check_solution_scipy,
)
from .test_min_cost_flow import (
    load_graph2_networkx,
    load_graph2_pandas,
    load_graph2_scipy,
)


class TestMaxFlow(unittest.TestCase):
    def setUp(self):
        self.expected_max_flow = 3.0

    def test_pandas(self):
        edge_data, _ = datasets.simple_graph_pandas()
        obj, sol = max_flow(edge_data, 0, 5)
        sol = sol[sol > 0]
        self.assertEqual(obj, self.expected_max_flow)
        candidate = {
            (0, 1): 1.0,
            (0, 2): 2.0,
            (1, 3): 1.0,
            (2, 3): 1.0,
            (2, 4): 1.0,
            (3, 5): 2.0,
            (4, 5): 1.0,
        }
        candidate2 = {
            (0, 1): 1.0,
            (0, 2): 2.0,
            (1, 3): 1.0,
            (2, 4): 2.0,
            (3, 5): 1.0,
            (4, 5): 2.0,
        }
        self.assertTrue(check_solution_pandas(sol, [candidate, candidate2]))

    def test_empty_pandas(self):
        edge_data, _ = datasets.simple_graph_pandas()
        edge_data["capacity"] = [0] * len(edge_data)
        obj, sol = max_flow(edge_data, 0, 5)
        self.assertEqual(obj, 0.0)

    def test_scipy(self):
        G, capacity, _, _ = datasets.simple_graph_scipy()
        G.data = capacity.data
        obj, sol = max_flow(G, 0, 5)
        self.assertEqual(obj, self.expected_max_flow)
        expected = np.array(
            [
                [0, 1, 2, 0, 0, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 2, 0],
                [0, 0, 0, 0, 0, 1],
                [0, 0, 0, 0, 0, 2],
            ]
        )
        self.assertTrue(check_solution_scipy(sol, [expected]))

    def test_empty_scipy(self):
        G, capacity, _, _ = datasets.simple_graph_scipy()
        # Use vsmall values to set the capacity to nearly zero
        G.data = np.repeat(1e-20, len(G.data))
        obj, sol = max_flow(G, 0, 5)
        self.assertEqual(obj, 0.0)

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_networkx(self):
        G = datasets.simple_graph_networkx()
        obj, sol = max_flow(G, 0, 5)
        self.assertEqual(obj, self.expected_max_flow)
        candidate = {
            (0, 1): {"flow": 1},
            (0, 2): {"flow": 2},
            (1, 3): {"flow": 1},
            (2, 4): {"flow": 2},
            (3, 5): {"flow": 1},
            (4, 5): {"flow": 2},
        }
        candidate2 = {
            (0, 1): {"flow": 1.0},
            (0, 2): {"flow": 2.0},
            (1, 3): {"flow": 1.0},
            (2, 3): {"flow": 1.0},
            (2, 4): {"flow": 1.0},
            (3, 5): {"flow": 2.0},
            (4, 5): {"flow": 1.0},
        }
        self.assertTrue(check_solution_networkx(sol, [candidate, candidate2]))

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_empty_networkx(self):
        G = datasets.simple_graph_networkx()
        nx.set_edge_attributes(G, 0, "capacity")
        obj, sol = max_flow(G, 0, 5)
        self.assertEqual(obj, 0.0)


class TestMaxFlow2(unittest.TestCase):
    def setUp(self):
        self.expected_max_flow = 23.0

    def test_pandas(self):
        edge_data, _ = load_graph2_pandas()
        obj, sol = max_flow(edge_data, 0, 4)
        sol = sol[sol > 0]
        self.assertEqual(obj, self.expected_max_flow)
        candidate = {
            (0, 1): 15.0,
            (0, 2): 8.0,
            (1, 3): 4.0,
            (1, 2): 1.0,
            (1, 4): 10.0,
            (2, 3): 4.0,
            (2, 4): 5.0,
            (3, 4): 8.0,
        }
        self.assertTrue(check_solution_pandas(sol, [candidate]))

    def test_scipy(self):
        G, capacity, _, _ = load_graph2_scipy()
        G.data = capacity.data
        obj, sol = max_flow(G, 0, 4)
        self.assertEqual(obj, self.expected_max_flow)
        expected = np.array(
            [
                [0, 15, 8, 0, 0],
                [0, 0, 1, 4, 10],
                [0, 0, 0, 9, 0],
                [0, 0, 0, 0, 13],
            ]
        )
        self.assertTrue(check_solution_scipy(sol, [expected]))

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_networkx(self):
        G = load_graph2_networkx()
        obj, sol = max_flow(G, 0, 4)
        self.assertEqual(obj, self.expected_max_flow)
        candidate = {
            (0, 1): {"flow": 15.0},
            (0, 2): {"flow": 8.0},
            (1, 3): {"flow": 4.0},
            (1, 2): {"flow": 1.0},
            (1, 4): {"flow": 10.0},
            (2, 3): {"flow": 4.0},
            (2, 4): {"flow": 5.0},
            (3, 4): {"flow": 8.0},
        }
        candidate2 = {
            (0, 1): {"flow": 15},
            (0, 2): {"flow": 8},
            (1, 2): {"flow": 1},
            (1, 3): {"flow": 4},
            (1, 4): {"flow": 10},
            (2, 3): {"flow": 9},
            (3, 4): {"flow": 13},
        }
        self.assertTrue(check_solution_networkx(sol, [candidate, candidate2]))
