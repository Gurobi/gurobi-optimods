import unittest

from gurobi_optimods.network_design import solve_network_design
from gurobi_optimods.datasets import load_network_design

try:
    import networkx as nx
except ImportError:
    nx = None


def create_tiny_graph():
    G = nx.DiGraph()
    G.add_edge(0, 1, capacity=1, fixed_cost=1, flow_cost=1)
    return G


@unittest.skipIf(nx is None, "networkx is not installed")
class TestNetworkDesign(unittest.TestCase):
    def test_empty(self):
        # Not sure what should happen for an empty graph
        G = nx.DiGraph()
        solve_network_design(G, {})

    def test_tiny_possible(self):
        # Tiny graph where the only answer is to build the edge (0, 1)
        G = create_tiny_graph()
        commodities = {0: {"Origin": 0, "Destination": 1, "Demand": 1}}
        sol, graph = solve_network_design(G, commodities)
        self.assertEqual(sol, 2)
        self.assertEqual(list(graph.edges(data=True)), [(0, 1, {"flow": {0: 1}})])

    def test_tiny_flow_exceeds_capacity(self):
        # Trivial graph but commodity quantity exceeds capacity
        G = create_tiny_graph()
        commodities = {0: {"Origin": 0, "Destination": 1, "Demand": 2}}
        with self.assertRaisesRegex(ValueError, "Unsatisfiable flows"):
            solve_network_design(G, commodities)

    def test_tiny_no_path_exists(self):
        # Trivial graph but no path exists for commodity from (1, 0)
        G = create_tiny_graph()
        commodities = {0: {"Origin": 1, "Destination": 0, "Demand": 1}}
        with self.assertRaisesRegex(ValueError, "Unsatisfiable flows"):
            solve_network_design(G, commodities)

    def test_origin_does_not_exist(self):
        G = create_tiny_graph()
        commodities = {0: {"Origin": 2, "Destination": 1, "Demand": 1}}
        with self.assertRaisesRegex(AssertionError, "Origin*"):
            solve_network_design(G, commodities)

    def test_destination_does_not_exist(self):
        G = create_tiny_graph()
        commodities = {0: {"Origin": 0, "Destination": 2, "Demand": 1}}
        with self.assertRaisesRegex(AssertionError, "Destination*"):
            solve_network_design(G, commodities)

    def test_network_flow(self):
        G = load_network_design()
        commodities = {
            0: {"Origin": 0, "Destination": 4, "Demand": 10},
            1: {"Origin": 2, "Destination": 4, "Demand": 15},
        }
        sol, graph = solve_network_design(G, commodities)
        self.assertEqual(sol, 176)
        print(list(graph.edges(data=True)))
        # self.assertEqual(
        #     list(graph.edges(data=True)),
        #     [
        #         (0, 1, {"flow": {0: 1.0, 1: 0.0}}),
        #         (1, 4, {"flow": {0: 1.0, 1: 0.0}}),
        #         (2, 3, {"flow": {0: 0.0, 1: 2.0 / 3}}),
        #         (2, 4, {"flow": {0: 0.0, 1: 1.0 / 3}}),
        #         (3, 4, {"flow": {0: 0.0, 1: 2.0 / 3}}),
        #     ],
        # )
