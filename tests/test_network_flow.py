# Unit tests for your network_flow. This is separate from the example tests, it should
# test smaller pieces of functionality, and can also test multiple cases.
# This template should be copied to tests/test_<network_flow-name>.py

import unittest

import scipy.sparse as sp

from gurobi_optimods.datasets import load_network_flow_example_data

from gurobi_optimods.network_flow import (
    network_flow,
    shortest_path,
    min_cost_flow,
    max_flow,
    min_cut,
)


class TestNetworkFlow(unittest.TestCase):
    def test_min_cost_flow(self):
        G, _, _ = load_network_flow_example_data()
        sol, cost = min_cost_flow(G)
        print(sol)
        self.assertEqual(cost, 150)
        self.assertEqual(
            [s[0] for s in sol],
            [
                (0, 1),
                (0, 2),
                (1, 3),
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 4),
            ],
        )
        self.assertEqual(
            [s[1] for s in sol],
            [12.0, 8.0, 4.0, 8.0, 11.0, 5.0, 10.0],
        )

    def test_shortest_path(self):
        G, source, sink = load_network_flow_example_data()
        sol, cost = shortest_path(G, source, sink)
        self.assertEqual(cost, 7)
        self.assertEqual([s[0] for s in sol], [(0, 2), (2, 3), (3, 4)])
        self.assertEqual([s[1] for s in sol], [1.0] * 3)

    def test_max_flow(self):
        G, source, sink = load_network_flow_example_data()
        sol, maxflow = max_flow(G, source, sink)
        self.assertEqual(maxflow, -23.0)
        self.assertEqual(
            [s[0] for s in sol],
            [(0, 1), (0, 2), (1, 3), (1, 2), (1, 4), (2, 3), (3, 4)],
        )
        self.assertEqual([s[1] for s in sol], [15.0, 8.0, 4.0, 1.0, 10.0, 9.0, 13.0])
