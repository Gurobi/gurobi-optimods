import io
import unittest

import pandas as pd
import scipy.sparse as sp

import gurobi_optimods.shortest_path as sp

node_data_csv = """
i,cost,time
0,1,5
1,1,4
2,1,3
3,1,2
4,1,1
5,1,0
"""

arc_data_csv = """
i,j,cost
0,1,1
0,2,3
1,2,-5
1,3,6
2,3,4
2,4,2
3,4,-5
3,5,1
4,2,0
4,5,3
"""


def load_graph_pandas():
    return (
        pd.read_csv(io.StringIO(node_data_csv)).set_index(["i"]),
        pd.read_csv(io.StringIO(arc_data_csv)).set_index(["i", "j"]),
    )


class TestShortestPath(unittest.TestCase):
    def test_negative_costs(self):
        node_data, arc_data = load_graph_pandas()
        path = sp.shortest_path(node_data, arc_data, 0, 5)
        self.assertEqual(path.cost, 4)
        arc_list = [0, 1, 2, 3, 4, 5]
        self.assertIsInstance(path.nodes, list)
        self.assertEqual(path.nodes, arc_list)

    # def test_infeasible(self):
    #     edge_data, node_data = datasets.simple_graph_pandas()
    #     # Add a node requesting more flow than is available.
    #     node_data["demand"].values[-1] = 10.0
    #     with self.assertRaisesRegex(ValueError, "Unsatisfiable flows"):
    #         obj, sol = mcf.min_cost_flow_pandas(edge_data, node_data)
