import io
import unittest

import numpy as np
import pandas as pd
import scipy.sparse as sp

try:
    import networkx as nx
except ImportError:
    nx = None

import gurobi_optimods.datasets as datasets
import gurobi_optimods.min_cost_flow as mcf

from .test_graph_utils import (
    check_solution_networkx,
    check_solution_networkx_multi,
    check_solution_pandas,
    check_solution_pandas_multi,
    check_solution_scipy,
)

edge_data2 = """
source,target,capacity,cost
0,1,15,4
0,2,8,4
1,3,4,2
1,2,20,2
1,4,10,6
2,3,15,1
2,4,5,3
3,4,20,2
"""


node_data2 = """
,demand
0,-20
1,0
2,0
3,5
4,15
"""

edge_data3 = """
source,target,capacity,cost
0,1,15,4
0,2,8,4
1,3,4,2
1,2,20,2
1,4,10,6
2,3,15,1
2,4,5,3
3,4,20,2
2,3,1,-100
2,3,1,10
"""


def load_graph2_pandas():
    return (
        pd.read_csv(io.StringIO(edge_data2)).set_index(["source", "target"]),
        pd.read_csv(io.StringIO(node_data2), index_col=0),
    )


def load_graph2_networkx(digraph=nx.DiGraph):
    edge_data, node_data = load_graph2_pandas()
    return datasets._convert_pandas_to_digraph(edge_data, node_data, digraph=digraph)


def load_graph2_scipy():
    edge_data, node_data = load_graph2_pandas()
    return datasets._convert_pandas_to_scipy(edge_data, node_data)


def load_graph3_pandas():
    return (
        pd.read_csv(io.StringIO(edge_data3)).set_index(["source", "target"]),
        pd.read_csv(io.StringIO(node_data2), index_col=0),
    )


def load_graph3_networkx(digraph=nx.DiGraph):
    edge_data, node_data = load_graph3_pandas()
    return datasets._convert_pandas_to_digraph(edge_data, node_data, digraph=digraph)


class TestMinCostFlow(unittest.TestCase):
    def test_pandas(self):
        edge_data, node_data = datasets.simple_graph_pandas()
        cost, sol = mcf.min_cost_flow_pandas(edge_data, node_data)
        sol = sol[sol > 0]
        self.assertEqual(cost, 31)
        candidate = {(0, 1): 1.0, (0, 2): 1.0, (1, 3): 1.0, (2, 4): 2.0, (4, 5): 2.0}
        self.assertIsInstance(sol, pd.Series)
        self.assertTrue(check_solution_pandas(sol, [candidate]))

    def test_infeasible(self):
        edge_data, node_data = datasets.simple_graph_pandas()
        # Add a node requesting more flow than is available.
        node_data["demand"].values[-1] = 10.0
        with self.assertRaisesRegex(ValueError, "Unsatisfiable flows"):
            obj, sol = mcf.min_cost_flow_pandas(edge_data, node_data)

    def test_scipy(self):
        G, cap, cost, demands = datasets.simple_graph_scipy()
        cost, sol = mcf.min_cost_flow_scipy(G, cap, cost, demands)
        self.assertEqual(cost, 31)
        candidate = np.array(
            [
                [0.0, 1.0, 1.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 2.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 2.0],
            ]
        )
        self.assertTrue(sp.issparse(sol))
        self.assertTrue(check_solution_scipy(sol, [candidate]))

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_networkx(self):
        G = datasets.simple_graph_networkx()
        cost, sol = mcf.min_cost_flow_networkx(G)
        self.assertEqual(cost, 31)
        expected = {
            (0, 1): {"flow": 1.0},
            (0, 2): {"flow": 1.0},
            (1, 3): {"flow": 1.0},
            (2, 4): {"flow": 2.0},
            (4, 5): {"flow": 2.0},
        }
        self.assertIsInstance(sol, nx.Graph)
        self.assertTrue(check_solution_networkx(sol, [expected]))

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_networkx_renamed(self):
        G = datasets.simple_graph_networkx()
        G = nx.relabel_nodes(G, {0: "s", 5: "t"})
        cost, sol = mcf.min_cost_flow_networkx(G)
        self.assertEqual(cost, 31)
        expected = {
            ("s", 1): {"flow": 1.0},
            ("s", 2): {"flow": 1.0},
            (1, 3): {"flow": 1.0},
            (2, 4): {"flow": 2.0},
            (4, "t"): {"flow": 2.0},
        }
        self.assertIsInstance(sol, nx.Graph)
        self.assertTrue(check_solution_networkx(sol, [expected]))


class TestMinCostFlow2(unittest.TestCase):
    def test_pandas(self):
        edge_data, node_data = load_graph2_pandas()
        cost, sol = mcf.min_cost_flow_pandas(edge_data, node_data)
        sol = sol[sol > 0]
        self.assertEqual(cost, 150)
        candidate = {
            (0, 1): 12.0,
            (0, 2): 8.0,
            (1, 3): 4.0,
            (1, 2): 8.0,
            (2, 3): 15.0,
            (2, 4): 1.0,
            (3, 4): 14.0,
        }
        candidate2 = {
            (0, 1): 12.0,
            (0, 2): 8.0,
            (1, 3): 4.0,
            (1, 2): 8.0,
            (2, 3): 11.0,
            (2, 4): 5.0,
            (3, 4): 10.0,
        }
        self.assertTrue(check_solution_pandas(sol, [candidate, candidate2]))

    def test_scipy(self):
        G, cap, cost, demands = load_graph2_scipy()
        cost, sol = mcf.min_cost_flow_scipy(G, cap, cost, demands)
        self.assertEqual(cost, 150)
        expected = np.array(
            [
                [0.0, 12.0, 8.0, 0.0, 0.0],
                [0.0, 0.0, 8.0, 4.0, 0.0],
                [0.0, 0.0, 0.0, 11.0, 5.0],
                [0.0, 0.0, 0.0, 0.0, 10.0],
            ]
        )
        self.assertTrue(check_solution_scipy(sol, [expected]))

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_networkx(self):
        G = load_graph2_networkx()
        cost, sol = mcf.min_cost_flow_networkx(G)
        self.assertEqual(cost, 150)
        candidate = {
            (0, 1): {"flow": 12.0},
            (0, 2): {"flow": 8.0},
            (1, 2): {"flow": 8.0},
            (1, 3): {"flow": 4.0},
            (2, 3): {"flow": 11.0},
            (2, 4): {"flow": 5.0},
            (3, 4): {"flow": 10.0},
        }
        candidate2 = {
            (0, 1): {"flow": 12.0},
            (0, 2): {"flow": 8.0},
            (1, 3): {"flow": 4.0},
            (1, 2): {"flow": 8.0},
            (2, 3): {"flow": 15.0},
            (2, 4): {"flow": 1.0},
            (3, 4): {"flow": 14.0},
        }
        self.assertTrue(check_solution_networkx(sol, [candidate, candidate2]))

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_networkx_renamed(self):
        G = load_graph2_networkx()
        G = nx.relabel_nodes(G, {0: "s", 4: "t"})
        cost, sol = mcf.min_cost_flow_networkx(G)
        self.assertEqual(cost, 150)
        candidate = {
            ("s", 1): {"flow": 12.0},
            ("s", 2): {"flow": 8.0},
            (1, 2): {"flow": 8.0},
            (1, 3): {"flow": 4.0},
            (2, 3): {"flow": 11.0},
            (2, "t"): {"flow": 5.0},
            (3, "t"): {"flow": 10.0},
        }
        candidate2 = {
            ("s", 1): {"flow": 12.0},
            ("s", 2): {"flow": 8.0},
            (1, 3): {"flow": 4.0},
            (1, 2): {"flow": 8.0},
            (2, 3): {"flow": 15.0},
            (2, "t"): {"flow": 1.0},
            (3, "t"): {"flow": 14.0},
        }
        self.assertTrue(check_solution_networkx(sol, [candidate, candidate2]))


class TestMinCostFlow3(unittest.TestCase):
    def test_pandas(self):
        edge_data, node_data = load_graph3_pandas()
        cost, sol = mcf.min_cost_flow_pandas(edge_data, node_data)
        sol = sol[sol > 0]
        self.assertEqual(cost, 49.0)

        candidate = pd.DataFrame(
            {
                "source": [0, 0, 1, 1, 2, 2, 3, 2],
                "target": [1, 2, 3, 2, 3, 4, 4, 3],
                "flow": [12.0, 8.0, 4.0, 8.0, 10.0, 5.0, 10.0, 1.0],
            }
        ).set_index(["source", "target"])

        self.assertTrue(check_solution_pandas_multi(sol, [candidate]))

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_networkx(self):
        G = load_graph3_networkx(digraph=nx.MultiDiGraph)
        cost, sol = mcf.min_cost_flow_networkx(G)
        self.assertEqual(cost, 49.0)
        candidate = [
            (0, 1, {"flow": 12.0}),
            (0, 2, {"flow": 8.0}),
            (1, 3, {"flow": 4.0}),
            (1, 2, {"flow": 8.0}),
            (2, 3, {"flow": 10.0}),
            (2, 3, {"flow": 1.0}),
            (2, 4, {"flow": 5.0}),
            (3, 4, {"flow": 10.0}),
        ]

        self.assertTrue(check_solution_networkx_multi(sol, [candidate]))
