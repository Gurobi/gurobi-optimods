import io
import unittest

import pandas as pd

try:
    import networkx as nx
except ImportError:
    nx = None

from gurobi_optimods.metromap import metromap

# to run only the unittests for the metromap optimod use this command
# python -m unittest tests.test_metromap.Testmap

small1_nodes = """
number
1
2
3
4
5
6
"""

small1_edges = """
source,target
1,2
2,5
3,4
2,6
5,3
6,3
2,3
"""
small1_linepath = """
linename,edge_source,edge_target
A,1,2
A,2,5
A,5,3
A,3,4
B,1,2
B,2,3
B,3,4
C,1,2
C,2,6
C,6,3
C,3,4
"""

small2_nodes = """
number,posx,posy
1,0,3
2,0,2
3,0,1
4,1,1
5,1,3
"""

small2_edges = """
source,target
1,2
2,1
1,5
5,1
2,3
3,2
3,4
4,3
4,5
5,4
"""


@unittest.skipIf(nx is None, "networkx is not installed")
class Testmap(unittest.TestCase):
    def test_node_degree_to_high(self):
        # create a graph one node in the middle and 9 adjacent nodes
        # node degree is 9
        G = nx.Graph()
        G.add_edge(0, 1)
        G.add_edge(0, 2)
        G.add_edge(0, 3)
        G.add_edge(0, 4)
        G.add_edge(0, 5)
        G.add_edge(0, 6)
        G.add_edge(0, 7)
        G.add_edge(0, 8)
        G.add_edge(0, 9)

        with self.assertRaises(ValueError) as error:
            graph_out, edge_directions = metromap(G)

        self.assertEqual(
            str(error.exception),
            "Node with number 0 has node degree 9. "
            "Octilinear representation is not possible for node degree larger than 8",
        )

    def test_small_no_geodata_no_lines(self):
        edge_data = pd.read_csv(io.StringIO(small1_edges))
        node_data = pd.read_csv(io.StringIO(small1_nodes))
        # create graph
        graph = nx.from_pandas_edgelist(
            edge_data.reset_index(), create_using=nx.Graph()
        )
        graph_out, edge_directions = metromap(graph)
        # a solution was computed, i.e., the graph is not empty
        self.assertEqual(graph_out.number_of_nodes(), 6)
        self.assertEqual(graph_out.number_of_edges(), 7)
        pos = nx.get_node_attributes(graph_out, "pos_oct")
        # also the attribute of the octilinear node positions is filled
        self.assertEqual(len(pos), 6)

    def test_small_no_geodata_with_lines(self):
        edge_data = pd.read_csv(io.StringIO(small1_edges))
        node_data = pd.read_csv(io.StringIO(small1_nodes))
        # create graph
        graph = nx.from_pandas_edgelist(
            edge_data.reset_index(), create_using=nx.Graph()
        )
        linepath_data = pd.read_csv(io.StringIO(small1_linepath))
        graph_out, edge_directions = metromap(graph, linepath_data)
        # a solution was computed, i.e., the graph is not empty
        self.assertEqual(graph_out.number_of_nodes(), 6)
        self.assertEqual(graph_out.number_of_edges(), 7)
        pos = nx.get_node_attributes(graph_out, "pos_oct")
        # also the attribute of the octilinear node positions is filled
        self.assertEqual(len(pos), 6)

    def test_small_geodata_no_lines(self):
        edge_data = pd.read_csv(io.StringIO(small2_edges))
        node_data = pd.read_csv(io.StringIO(small2_nodes))
        # create graph
        graph = nx.from_pandas_edgelist(
            edge_data.reset_index(), create_using=nx.Graph()
        )
        # add x-, y-coordinates as node attribute
        for number, row in node_data.set_index("number").iterrows():
            graph.add_node(number, pos=(row["posx"], row["posy"]))

        graph_out, edge_directions = metromap(graph)
        # a solution was computed, i.e., the graph is not empty
        self.assertEqual(graph_out.number_of_nodes(), 5)
        self.assertEqual(graph_out.number_of_edges(), 5)
        pos = nx.get_node_attributes(graph_out, "pos_oct")
        # also the attribute of the octilinear node positions is filled
        self.assertEqual(len(pos), 5)

    def test_small_geodata_min_dist(self):
        edge_data = pd.read_csv(io.StringIO(small2_edges))
        node_data = pd.read_csv(io.StringIO(small2_nodes))
        # create graph
        graph = nx.from_pandas_edgelist(
            edge_data.reset_index(), create_using=nx.Graph()
        )
        # add x-, y-coordinates as node attribute
        for number, row in node_data.set_index("number").iterrows():
            graph.add_node(number, pos=(row["posx"], row["posy"]))

        graph_out, edge_directions = metromap(
            graph, penalty_edge_directions=0, penalty_distance=1
        )
        # a solution was computed, i.e., the graph is not empty
        self.assertEqual(graph_out.number_of_nodes(), 5)
        self.assertEqual(graph_out.number_of_edges(), 5)
        pos = nx.get_node_attributes(graph_out, "pos_oct")
        self.assertTrue(edge_directions[(3, 4)] == 1 or edge_directions[(1, 5)] == 7)

    def test_small_geodata_min_geodiff(self):
        edge_data = pd.read_csv(io.StringIO(small2_edges))
        node_data = pd.read_csv(io.StringIO(small2_nodes))
        # create graph
        graph = nx.from_pandas_edgelist(
            edge_data.reset_index(), create_using=nx.Graph()
        )
        # add x-, y-coordinates as node attribute
        for number, row in node_data.set_index("number").iterrows():
            graph.add_node(number, pos=(row["posx"], row["posy"]))

        graph_out, edge_directions = metromap(
            graph, penalty_edge_directions=1, penalty_distance=0
        )
        # a solution was computed, i.e., the graph is not empty
        self.assertEqual(graph_out.number_of_nodes(), 5)
        self.assertEqual(graph_out.number_of_edges(), 5)
        pos = nx.get_node_attributes(graph_out, "pos_oct")
        self.assertTrue(edge_directions[(3, 4)] == 0)
        self.assertTrue(edge_directions[(1, 5)] == 0)
