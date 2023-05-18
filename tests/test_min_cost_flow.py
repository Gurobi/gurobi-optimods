import unittest

import gurobipy as gp
import numpy as np
from gurobipy import GRB

try:
    import networkx as nx
except ImportError:
    nx = None

from numpy.testing import assert_array_equal, assert_allclose

import gurobi_optimods.datasets as datasets
import gurobi_optimods.min_cost_flow as mcf
from gurobi_optimods.min_cost_flow import solve_min_cost_flow


class TestMinCostFlow(unittest.TestCase):
    def test_pandas(self):
        edge_data, node_data = datasets.load_graph()
        cost, sol = mcf.min_cost_flow(edge_data, node_data)
        sol = sol[sol > 0]
        self.assertEqual(cost, 31)
        self.assertEqual(sol.tolist(), [1.0, 1.0, 1.0, 2.0, 2.0])
        self.assertEqual(sol.index.tolist(), [(0, 1), (0, 2), (1, 3), (2, 4), (4, 5)])

    def test_infeasible(self):
        edge_data, node_data = datasets.load_graph()
        # Add a node requesting more flow than is available.
        node_data["demand"].values[-1] = 10.0
        with self.assertRaisesRegex(ValueError, "Unsatisfiable flows"):
            obj, sol = mcf.min_cost_flow(edge_data, node_data)

    def test_scipy(self):
        G, cap, cost, demands = datasets.load_graph_scipy()
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
    def test_networkx(self):
        G = datasets.load_graph_networkx()
        cost, sol = mcf.min_cost_flow_networkx(G)
        self.assertEqual(cost, 31)
        self.assertEqual(
            sol, {(0, 1): 1.0, (0, 2): 1.0, (1, 3): 1.0, (2, 4): 2.0, (4, 5): 2.0}
        )


class TestMinCostFlow2(unittest.TestCase):
    def test_pandas(self):
        edge_data, node_data = datasets.load_graph2()
        cost, sol = mcf.min_cost_flow(edge_data, node_data)
        sol = sol[sol > 0]
        self.assertEqual(cost, 150)
        self.assertEqual(sol.tolist(), [12.0, 8.0, 4.0, 8.0, 11.0, 5.0, 10.0])
        self.assertEqual(
            sol.index.tolist(),
            [(0, 1), (0, 2), (1, 3), (1, 2), (2, 3), (2, 4), (3, 4)],
        )

    def test_scipy(self):
        G, cap, cost, demands = datasets.load_graph2_scipy()
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
    def test_networkx(self):
        G = datasets.load_graph2_networkx()
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


class TestSolveMinCostFlowUtil(unittest.TestCase):
    # FIXME repurpose these to target the mod

    def setUp(self):
        self.env = gp.Env()

    def tearDown(self):
        self.env.close()

    def test_max_flow(self):
        # Max flow network from a small bipartite matching model, augmented
        # with a sink->source edge to create a min-cost flow
        cost, flows = solve_min_cost_flow(
            env=self.env,
            edge_source=np.array([5, 5, 0, 1, 1, 2, 3, 4, 6]),
            edge_target=np.array([0, 1, 4, 2, 4, 6, 6, 6, 5]),
            capacity=np.array([1.0] * 8 + [GRB.INFINITY]),
            cost=np.array([0.0] * 8 + [-1.0]),
            demand=np.zeros(7),
        )
        # Max flow is 2. Any optimal solution has 2 units along the sink->source
        # arc and 2 units in each of the 3 network layers.
        self.assertEqual(cost, -2.0)
        self.assertEqual(flows[-1], 2.0)
        self.assertEqual(flows.sum(), 8.0)

    def test_linear(self):
        # Silly example to get the demand direction right
        cost, flows = solve_min_cost_flow(
            env=self.env,
            edge_source=np.array([0, 1, 2, 3]),
            edge_target=np.array([1, 2, 3, 4]),
            capacity=np.ones(4),
            cost=np.ones(4),
            demand=np.array([-1.0, 0.0, 0.0, 0.0, 1.0]),
        )
        self.assertEqual(cost, 4.0)
        assert_allclose(flows, np.ones(4))

    def test_infeasible(self):
        # A very silly example
        with self.assertRaisesRegex(ValueError, "Unsatisfiable flows"):
            solve_min_cost_flow(
                env=self.env,
                edge_source=np.array([0]),
                edge_target=np.array([1]),
                capacity=np.array([1.0]),
                cost=np.array([1.0]),
                demand=np.array([-1.0, 2.0]),
            )
