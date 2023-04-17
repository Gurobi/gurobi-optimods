import numpy as np

import gurobipy as gp
from gurobipy import GRB


def maximum_weighted_independent_set(adjacency_matrix, weights):
    """Find a set of mutually non-adjacent vertices with maximum weighted sum.

    :param adjacency_matrix: The upper triangular adjacency matrix.
    :type adjacency_matrix: :class:`sp.spmatrix`
    :param weights: Vertex weight array.
    :type weights: :class:`np.array`
    :return: The maximum weighted independent set array.
    :rtype: :class:`np.array`
    """
    with gp.Env() as env, gp.Model("mwis", env=env) as model:
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
