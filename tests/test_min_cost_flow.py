import unittest

import numpy as np
import pandas as pd
import scipy.sparse as sp

from numpy.testing import assert_array_equal

try:
    import networkx as nx
except ImportError:
    nx = None

from gurobi_optimods.datasets import load_min_cost_flow, load_min_cost_flow2
from gurobi_optimods.min_cost_flow import (
    min_cost_flow,
    min_cost_flow_networkx,
    min_cost_flow_scipy,
)


def _convert_pandas_to_digraph(edge_data, node_data):
    """
    Convert from a pandas DataFrame to a networkx.DiGraph with the appropriate
    attributes. For edges: `capacity`, and `cost`. For nodes: `demand`.
    """
    G = nx.from_pandas_edgelist(
        edge_data.reset_index(), create_using=nx.DiGraph(), edge_attr=True
    )
    for i, d in node_data.iterrows():
        G.add_node(i, demand=d.demand)
    return G


def _convert_pandas_to_scipy(edge_data, node_data):
    """
    Convert from a pandas DataFrame to several scipy.sparse.coo_matrix contain
    the graph structure, the capacity and cost values per edge, and the demand
    values per node.
    """
    coords = edge_data.index.to_numpy()

    a0 = np.array([c[0] for c in coords])
    a1 = np.array([c[1] for c in coords])

    data = np.ones(len(coords))
    G = sp.coo_matrix((data, (a0, a1)))

    data = edge_data["capacity"].values
    cap = sp.coo_matrix((data, (a0, a1)))

    data = edge_data["cost"].values
    costs = sp.coo_matrix((data, (a0, a1)))

    return G, cap, costs, node_data["demand"].values


class TestMinCostFlow(unittest.TestCase):
    def setUp(self):
        self.edge_data, self.node_data = load_min_cost_flow()

    def test_min_cost_flow(self):
        cost, sol = min_cost_flow(self.edge_data, self.node_data)
        self.assertEqual(cost, 150)
        self.assertEqual(sol.tolist(), [12.0, 8.0, 4.0, 8.0, 0.0, 11.0, 5.0, 10.0, 0.0])
        self.assertEqual(
            sol.index.tolist(),
            [(0, 1), (0, 2), (1, 3), (1, 2), (1, 4), (2, 3), (2, 4), (3, 4), (4, 2)],
        )

    def test_min_cost_flow_scipy(self):
        G, cap, cost, demands = _convert_pandas_to_scipy(self.edge_data, self.node_data)
        cost, sol = min_cost_flow_scipy(G, cap, cost, demands)
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
        G = _convert_pandas_to_digraph(self.edge_data, self.node_data)
        cost, sol = min_cost_flow_networkx(G)
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


class TestMinCostFlow2(unittest.TestCase):
    def setUp(self):
        self.edge_data, self.node_data = load_min_cost_flow2()

    def test_min_cost_flow(self):
        cost, sol = min_cost_flow(self.edge_data, self.node_data)
        self.assertEqual(cost, 31)
        self.assertEqual(sol.tolist(), [1.0, 1.0, 1.0, 0.0, 2.0, 0.0, 2.0])
        self.assertEqual(
            sol.index.tolist(), [(0, 1), (0, 2), (1, 3), (2, 3), (2, 4), (3, 5), (4, 5)]
        )

    def test_min_cost_flow_scipy(self):
        G, cap, cost, demands = _convert_pandas_to_scipy(self.edge_data, self.node_data)
        cost, sol = min_cost_flow_scipy(G, cap, cost, demands)
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
        G = _convert_pandas_to_digraph(self.edge_data, self.node_data)
        cost, sol = min_cost_flow_networkx(G)
        self.assertEqual(cost, 31)
        self.assertEqual(
            sol, {(0, 1): 1.0, (0, 2): 1.0, (1, 3): 1.0, (2, 4): 2.0, (4, 5): 2.0}
        )
