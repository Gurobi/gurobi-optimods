import unittest

import numpy as np

try:
    import networkx as nx
except ImportError:
    nx = None

import gurobi_optimods.datasets as datasets
from gurobi_optimods.max_flow import max_flow

from .test_graph_utils import (
    check_solution_pandas,
    check_solution_scipy,
    check_solution_networkx,
)


class TestMaxFlow(unittest.TestCase):
    def setUp(self):
        self.expected_max_flow = 3.0

    def test_infeasible(self):
        edge_data, _ = datasets.load_graph()
        edge_data["capacity"][(0, 1)] = -1.0
        with self.assertRaisesRegex(ValueError, "Unsatisfiable flows"):
            obj, sol = max_flow(edge_data, 0, 5)

    def test_pandas(self):
        edge_data, _ = datasets.load_graph()
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

    def test_scipy(self):
        G, capacity, _, _ = datasets.load_graph_scipy()
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

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_networkx(self):
        G = datasets.load_graph_networkx()
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


class TestMaxFlow2(unittest.TestCase):
    def setUp(self):
        self.expected_max_flow = 23.0

    def test_pandas(self):
        edge_data, _ = datasets.load_graph2()
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
        G, capacity, _, _ = datasets.load_graph2_scipy()
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
        G = datasets.load_graph2_networkx()
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
