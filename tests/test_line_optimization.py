import io
import pathlib
import unittest

import pandas as pd

import gurobi_optimods.datasets as datasets
import gurobi_optimods.line_optimization as lop

DATA_FILE_DIR = pathlib.Path(__file__).parent / "data"

edge_data2A = """
sourc,target,length,time
0,1,15,40
0,2,8,40
"""

edge_data2 = """
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

node_data2 = """
number
0
1
2
3
4
"""

lines2 = """
linename,capacity,fixCost,operatingCost
L1,20,9,4
L2,20,9,2
L3,20,9,3
"""

line_path2 = """
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

demand2 = """
source,target,demand
0,3,20
0,5,30
3,4,30
"""

demand2A = """
source,target,demand
0,4,-30
0,3,10
2,1,10
"""

# python -m unittest tests.test_line_optimization.Testlop


class Testlop(unittest.TestCase):
    def test_sol(self):
        (
            node_data,
            edge_data,
            line_data,
            linepath_data,
            demand_data,
        ) = datasets.load_siouxfalls_network_data()
        frequencies = [1, 3]
        objCost, finalLines = lop.line_optimization(
            node_data,
            edge_data,
            line_data,
            linepath_data,
            demand_data,
            frequencies,
            True,
        )
        self.assertEqual(objCost, 211)
        self.assertEqual(len(finalLines), 8)

    def test_wrong_data_format(self):
        edge_data = pd.read_csv(io.StringIO(edge_data2A))
        node_data = pd.read_csv(io.StringIO(node_data2))
        linepath_data = pd.read_csv(io.StringIO(line_path2))
        line_data = pd.read_csv(io.StringIO(lines2))
        demand_data = pd.read_csv(io.StringIO(demand2))
        frequencies = [3, 6]
        objCost, finalLines = lop.line_optimization(
            node_data, edge_data, line_data, linepath_data, demand_data, frequencies
        )
        self.assertEqual(objCost, -1)
        self.assertEqual(len(finalLines), 0)

    def test_wrong_data_format2(self):
        edge_data = pd.read_csv(io.StringIO(edge_data2))
        node_data = pd.read_csv(io.StringIO(node_data2))
        linepath_data = pd.read_csv(io.StringIO(line_path2))
        line_data = pd.read_csv(io.StringIO(lines2))
        demand_data = pd.read_csv(io.StringIO(demand2A))
        frequencies = [3, 6, 9, 18]
        objCost, finalLines = lop.line_optimization(
            node_data, edge_data, line_data, linepath_data, demand_data, frequencies
        )
        self.assertEqual(objCost, -1)
        self.assertEqual(len(finalLines), 0)

    def test_wrong_data_format3(self):
        edge_data = pd.read_csv(io.StringIO(edge_data2))
        node_data = pd.read_csv(io.StringIO(node_data2))
        linepath_data = pd.read_csv(io.StringIO(line_path2))
        line_data = pd.read_csv(io.StringIO(lines2))
        demand_data = pd.read_csv(io.StringIO(demand2))
        frequencies = [1, -2, 3]
        objCost, finalLines = lop.line_optimization(
            node_data, edge_data, line_data, linepath_data, demand_data, frequencies
        )

    def test_smallAllPaths(self):
        edge_data = pd.read_csv(io.StringIO(edge_data2))
        node_data = pd.read_csv(io.StringIO(node_data2))
        linepath_data = pd.read_csv(io.StringIO(line_path2))
        line_data = pd.read_csv(io.StringIO(lines2))
        demand_data = pd.read_csv(io.StringIO(demand2))
        frequencies = [1, 2, 3]
        objCost, finalLines = lop.line_optimization(
            node_data,
            edge_data,
            line_data,
            linepath_data,
            demand_data,
            frequencies,
            False,
        )
        self.assertEqual(objCost, 21)
        self.assertEqual(finalLines, [("L1", 3)])

    def test_shortestPath(self):
        edge_data = pd.read_csv(io.StringIO(edge_data2))
        node_data = pd.read_csv(io.StringIO(node_data2))
        linepath_data = pd.read_csv(io.StringIO(line_path2))
        line_data = pd.read_csv(io.StringIO(lines2))
        demand_data = pd.read_csv(io.StringIO(demand2))
        frequencies = [1, 2, 3]
        objCost, finalLines = lop.line_optimization(
            node_data,
            edge_data,
            line_data,
            linepath_data,
            demand_data,
            frequencies,
            True,
        )
        self.assertEqual(objCost, 50)
        self.assertEqual(finalLines, [("L1", 2), ("L2", 3), ("L3", 3)])
