import contextlib
import io
import pathlib
import unittest

import pandas as pd

try:
    import networkx as nx
except ImportError:
    nx = None

import gurobi_optimods.datasets as datasets
import gurobi_optimods.metromap as map

DATA_FILE_DIR = pathlib.Path(__file__).parent / "data"


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
    def test_small_no_geodata_no_lines(self):
        edge_data = pd.read_csv(io.StringIO(small1_edges))
        node_data = pd.read_csv(io.StringIO(small1_nodes))
        # create graph
        graph = nx.from_pandas_edgelist(
            edge_data.reset_index(), create_using=nx.Graph()
        )
        graph_out, edge_directions = map.metromap(graph)
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
        graph_out, edge_directions = map.metromap(graph, linepath_data)
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

        graph_out, edge_directions = map.metromap(graph)
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

        graph_out, edge_directions = map.metromap(
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

        graph_out, edge_directions = map.metromap(
            graph, penalty_edge_directions=1, penalty_distance=0
        )
        # a solution was computed, i.e., the graph is not empty
        self.assertEqual(graph_out.number_of_nodes(), 5)
        self.assertEqual(graph_out.number_of_edges(), 5)
        pos = nx.get_node_attributes(graph_out, "pos_oct")
        self.assertTrue(edge_directions[(3, 4)] == 0)
        self.assertTrue(edge_directions[(1, 5)] == 0)

    # def test_real_example(self):
    #     (graph, linepath_data) = datasets.load_sberlin_graph_data()
    #     graph_out, edge_directions = map.metromap(graph, include_planarity=False)
    #     print("print solution")
    #     map.plot_network(graph_out, edge_directions, linepath_data)
