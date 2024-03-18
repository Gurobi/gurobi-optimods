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
import gurobi_optimods.line_optimization as lop

DATA_FILE_DIR = pathlib.Path(__file__).parent / "data"

edges_wrong = """
sourc,target,length,time
0,1,15,40
0,2,8,40
"""

edges = """
source,target,time
0,1,10
0,2,10
1,3,30
2,3,20
2,4,50
3,4,40
3,5,20
4,5,30
"""

nodes = """
number
0
1
2
3
4
5
"""

nodes_unknown_node = """
number
0
1
2
3
4
5
6
"""

lines = """
linename,capacity,fixCost,operatingCost
L1,20,9,4
L2,20,9,2
L3,20,9,3
"""

linepath = """
linename,edgeSource,edgeTarget
L1,0,1
L1,1,3
L1,3,4
L1,4,5
L2,2,3
L2,3,5
L3,0,2
L3,2,4
L3,4,5
"""

demand = """
source,target,demand
0,3,20
0,5,30
3,4,30
"""

demand_negative = """
source,target,demand
0,4,-30
0,3,10
2,1,10
"""

demand_unknown_node = """
source,target,demand
0,3,20
0,5,30
3,4,30
2,6,30
"""

demand_no_connection = """
source,target,demand
2,1,20
"""

# python -m unittest tests.test_line_optimization.Testlop


class Testlop(unittest.TestCase):
    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_sol(self):
        (
            node_data,
            edge_data,
            line_data,
            linepath_data,
            demand_data,
        ) = datasets.load_siouxfalls_network_data()
        frequencies = [1, 3]
        obj_cost, final_lines = lop.line_optimization(
            node_data,
            edge_data,
            line_data,
            linepath_data,
            demand_data,
            frequencies,
            True,
        )
        self.assertEqual(obj_cost, 211)
        self.assertEqual(len(final_lines), 8)

    def test_wrong_data_format(self):
        edge_data = pd.read_csv(io.StringIO(edges_wrong))
        node_data = pd.read_csv(io.StringIO(nodes))
        linepath_data = pd.read_csv(io.StringIO(linepath))
        line_data = pd.read_csv(io.StringIO(lines))
        demand_data = pd.read_csv(io.StringIO(demand))
        frequencies = [3, 6]

        captured_output = io.StringIO()
        try:
            with contextlib.redirect_stdout(captured_output):
                lop.line_optimization(
                    node_data,
                    edge_data,
                    line_data,
                    linepath_data,
                    demand_data,
                    frequencies,
                )
        except ValueError:
            self.assertEqual(
                "column source not present in edge_data\n", captured_output.getvalue()
            )

    def test_wrong_data_neg_demand(self):
        edge_data = pd.read_csv(io.StringIO(edges))
        node_data = pd.read_csv(io.StringIO(nodes))
        linepath_data = pd.read_csv(io.StringIO(linepath))
        line_data = pd.read_csv(io.StringIO(lines))
        demand_data = pd.read_csv(io.StringIO(demand_negative))
        frequencies = [3, 18]
        captured_output = io.StringIO()
        try:
            with contextlib.redirect_stdout(captured_output):
                lop.line_optimization(
                    node_data,
                    edge_data,
                    line_data,
                    linepath_data,
                    demand_data,
                    frequencies,
                )
        except ValueError:
            self.assertEqual(
                "All demand should be non-negative! Some value is negative in demand_data\n",
                captured_output.getvalue(),
            )

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_wrong_data_node_not_found(self):
        edge_data = pd.read_csv(io.StringIO(edges))
        node_data = pd.read_csv(io.StringIO(nodes_unknown_node))
        linepath_data = pd.read_csv(io.StringIO(linepath))
        line_data = pd.read_csv(io.StringIO(lines))
        demand_data = pd.read_csv(io.StringIO(demand_unknown_node))
        frequencies = [3, 6, 9, 18]
        with self.assertRaises(ValueError) as error:
            lop.line_optimization(
                node_data,
                edge_data,
                line_data,
                linepath_data,
                demand_data,
                frequencies,
                True,
            )
        self.assertEqual(str(error.exception), "demand node 6 not found in edges")

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_wrong_data_no_path_found(self):
        edge_data = pd.read_csv(io.StringIO(edges))
        node_data = pd.read_csv(io.StringIO(nodes))
        linepath_data = pd.read_csv(io.StringIO(linepath))
        line_data = pd.read_csv(io.StringIO(lines))
        demand_data = pd.read_csv(io.StringIO(demand_no_connection))
        frequencies = [18]
        with self.assertRaises(ValueError) as error:
            lop.line_optimization(
                node_data,
                edge_data,
                line_data,
                linepath_data,
                demand_data,
                frequencies,
                True,
            )
        self.assertEqual(
            str(error.exception), "no path found for connection from 2 to 1"
        )

    def test_smallAllPaths(self):
        edge_data = pd.read_csv(io.StringIO(edges))
        node_data = pd.read_csv(io.StringIO(nodes))
        linepath_data = pd.read_csv(io.StringIO(linepath))
        line_data = pd.read_csv(io.StringIO(lines))
        demand_data = pd.read_csv(io.StringIO(demand))
        frequencies = [1, 2, 3]
        obj_cost, final_lines = lop.line_optimization(
            node_data,
            edge_data,
            line_data,
            linepath_data,
            demand_data,
            frequencies,
            False,
        )
        self.assertEqual(obj_cost, 21)
        self.assertEqual(final_lines, [("L1", 3)])

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_shortestPath(self):
        edge_data = pd.read_csv(io.StringIO(edges))
        node_data = pd.read_csv(io.StringIO(nodes))
        linepath_data = pd.read_csv(io.StringIO(linepath))
        line_data = pd.read_csv(io.StringIO(lines))
        demand_data = pd.read_csv(io.StringIO(demand))
        frequencies = [1, 2, 3]
        obj_cost, final_lines = lop.line_optimization(
            node_data,
            edge_data,
            line_data,
            linepath_data,
            demand_data,
            frequencies,
            True,
        )
        self.assertEqual(obj_cost, 50)
        self.assertEqual(final_lines, [("L1", 2), ("L2", 3), ("L3", 3)])
