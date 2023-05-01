# Unit tests for your mod. This is separate from the example tests, it should
# test smaller pieces of functionality, and can also test multiple cases.
# This template should be copied to tests/test_<mod-name>.py

import unittest

import scipy.sparse as sp
import networkx as nx

from gurobi_optimods.network_design import solve_network_design
from gurobi_optimods.datasets import load_network_design, load_commodities


class TestNetworkDesign(unittest.TestCase):
    def setUp(self):
        self.G = load_network_design()
        self.commodities = load_commodities()

    def test_network_dataset(self):
        self.assertEqual(self.G.number_of_nodes(), 5)
        self.assertEqual(self.G.number_of_edges(), 9)

    def test_network_flow(self):
        sol, cost = solve_network_design(self.G, self.commodities)
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
