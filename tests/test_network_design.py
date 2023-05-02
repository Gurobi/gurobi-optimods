# Unit tests for your mod. This is separate from the example tests, it should
# test smaller pieces of functionality, and can also test multiple cases.
# This template should be copied to tests/test_<mod-name>.py

import unittest

from gurobi_optimods.network_design import solve_network_design
from gurobi_optimods.datasets import load_network_design, load_commodities


class TestNetworkDesign(unittest.TestCase):
    def setUp(self):
        self.G = load_network_design()
        self.commodities = {
            0: {"Origin": 0, "Destination": 4, "Demand": 10},
            1: {"Origin": 2, "Destination": 4, "Demand": 15},
        }

    def test_network_dataset(self):
        self.assertEqual(self.G.number_of_nodes(), 5)
        self.assertEqual(self.G.number_of_edges(), 9)

    def test_commodities(self):
        for k in self.commodities.keys():
            assert k in self.G.nodes()

    def test_network_flow(self):
        sol, graph = solve_network_design(self.G, self.commodities)
        self.assertEqual(sol, 176)
