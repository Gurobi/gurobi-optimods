# Unit tests for your network_flow. This is separate from the example tests, it should
# test smaller pieces of functionality, and can also test multiple cases.
# This template should be copied to tests/test_<network_flow-name>.py

import unittest

import scipy.sparse as sp
import networkx as nx

from gurobi_optimods.datasets import (
    load_network_flow_example_data,
    load_network_flow_example_data2,
)

from gurobi_optimods.network_flow import (
    network_flow,
    shortest_path,
    min_cost_flow,
    max_flow,
    min_cut,
)


class TestNetworkFlow(unittest.TestCase):
    def setUp(self):
        self.G, self.source, self.sink = load_network_flow_example_data()

    def test_network_flow(self):
        sol, cost = network_flow(self.G)
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

    def test_min_cost_flow(self):
        sol, cost = min_cost_flow(self.G)
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

    def test_shortest_path(self):
        sol, cost = shortest_path(self.G, self.source, self.sink)
        self.assertEqual(cost, 7)
        # Multiple solutions
        self.assertIn(
            sol,
            [
                {(0, 2): 1, (2, 4): 1},
                {(0, 2): 1, (2, 3): 1, (3, 4): 1},
            ],
        )

    def test_max_flow(self):
        sol, maxflow = max_flow(self.G, self.source, self.sink)
        self.assertEqual(maxflow, 23.0)
        # Multiple solutions
        self.assertIn(
            sol,
            [
                {
                    (0, 1): 15.0,
                    (0, 2): 8.0,
                    (1, 3): 4.0,
                    (1, 2): 1.0,
                    (1, 4): 10.0,
                    (2, 3): 9.0,
                    (3, 4): 13.0,
                },
                {
                    (0, 1): 15.0,
                    (0, 2): 8.0,
                    (1, 3): 4.0,
                    (1, 2): 1.0,
                    (1, 4): 10.0,
                    (2, 3): 4.0,
                    (2, 4): 5.0,
                    (3, 4): 8.0,
                },
            ],
        )

    def test_min_cut(self):
        partition, flow = min_cut(self.G, self.source, self.sink)
        self.assertEqual(partition[0], set({0}))
        self.assertEqual(partition[1], set({1, 2, 3, 4}))

    def test_min_cut2(self):
        G = nx.DiGraph()
        G.add_edge("x", "a", capacity=3.0)
        G.add_edge("x", "b", capacity=1.0)
        G.add_edge("a", "c", capacity=3.0)
        G.add_edge("b", "c", capacity=5.0)
        G.add_edge("b", "d", capacity=4.0)
        G.add_edge("d", "e", capacity=2.0)
        G.add_edge("c", "y", capacity=2.0)
        G.add_edge("e", "y", capacity=3.0)
        partition, flow = min_cut(G, "x", "y")
        print(partition)
        self.assertEqual(partition[0], set({"c", "a", "x"}))
        self.assertEqual(partition[1], set({"d", "y", "e", "b"}))


class TestNetworkFlow2(unittest.TestCase):
    def setUp(self):
        self.G, self.source, self.sink = load_network_flow_example_data2()

    def test_network_flow(self):
        sol, cost = network_flow(self.G)
        self.assertEqual(
            sol, {("s", 1): 1.0, ("s", 2): 1.0, (1, 3): 1.0, (2, 4): 2.0, (4, "t"): 2.0}
        )
        self.assertEqual(cost, 31)

    def test_min_cost_flow(self):
        sol, cost = min_cost_flow(self.G)
        self.assertEqual(
            sol, {("s", 1): 1.0, ("s", 2): 1.0, (1, 3): 1.0, (2, 4): 2.0, (4, "t"): 2.0}
        )
        self.assertEqual(cost, 31)

    def test_max_flow(self):
        sol, flow = max_flow(self.G, self.source, self.sink)
        self.assertEqual(
            sol,
            {
                ("s", 1): 1.0,
                ("s", 2): 2.0,
                (1, 3): 1.0,
                (2, 4): 2.0,
                (3, "t"): 1.0,
                (4, "t"): 2.0,
            },
        )
        self.assertEqual(flow, 3.0)

    def test_min_cut(self):
        partition, flow = min_cut(self.G, self.source, self.sink)
        self.assertEqual(partition[0], set({"s", 1}))
        self.assertEqual(partition[1], set({2, 3, 4, "t"}))

    def test_shortest_path(self):
        sol, cost = shortest_path(self.G, self.source, self.sink)
        self.assertEqual(sol, {("s", 1): 1.0, (1, 3): 1.0, (3, "t"): 1.0})
        self.assertEqual(cost, 11)
