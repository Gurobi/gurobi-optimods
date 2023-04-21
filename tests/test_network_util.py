import unittest

import numpy as np
import gurobipy as gp
from gurobipy import GRB
from numpy.testing import assert_allclose

from gurobi_optimods.network_util import solve_min_cost_flow


class GurobiTestBase(unittest.TestCase):
    def setUp(self):
        self.env = gp.Env()

    def tearDown(self):
        self.env.close()


class TestSolveMinCostFlow(GurobiTestBase):
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
