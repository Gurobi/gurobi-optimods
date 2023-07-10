"""
Maximum Weighted Independent Set
--------------------------------
"""

import gurobipy as gp
import numpy as np
from gurobipy import GRB

from gurobi_optimods.utils import optimod


@optimod()
def maximum_weighted_independent_set(adjacency_matrix, weights, *, create_env):
    """Find a set of mutually non-adjacent vertices with maximum weighted sum.

    Parameters
    ----------
    adjacency_matrix : spmatrix
        The upper triangular adjacency matrix.
    weights : ndarray
        Vertex weight array.

    Returns
    -------
    ndarray
        The maximum weighted independent set array.
    """
    with create_env() as env, gp.Model("mwis", env=env) as model:
        # x_i: 1 if vertex i is in the independent set and 0 otherwise
        x = model.addMVar(len(weights), vtype=GRB.BINARY, name="x")
        # Maximize the sum of the vertex weights in the independent set
        model.setObjective(weights @ x, sense=GRB.MAXIMIZE)
        # The independent set contains non-adjacent vertices
        model.addConstrs(
            (
                x[i] + x[j] <= 1
                for i, j in zip(
                    adjacency_matrix.tocoo().row, adjacency_matrix.tocoo().col
                )
            ),
            name="no_adjacent_vertices",
        )
        model.optimize()
        return np.where(x.X >= 0.5)[0]
