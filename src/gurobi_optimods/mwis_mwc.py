"""
Maximum Weighted Independent Set/Clique
---------------------------------------
"""

from dataclasses import dataclass

import gurobipy as gp
import numpy as np
import scipy.sparse as sp
from gurobipy import GRB

from gurobi_optimods.utils import optimod


@dataclass
class Result:
    """
    Data class representing the maximum weighted independent set
    (clique) and its weight

    Attributes
    ----------
    x : ndarray
        The maximum weighted independent set (clique) array
    f : float
        The total weight of the maximum weighted independent set (clique)
    """

    x: np.ndarray
    f: float


def check_input(adjacency_matrix, weights):
    if not isinstance(adjacency_matrix, sp.spmatrix):
        raise ValueError(
            "The adjacency matrix should be a scipy sparse matrix in CSR format."
        )
    if not isinstance(weights, np.ndarray):
        raise ValueError("The weights of the vertices should be a numpy array.")

    if not np.allclose(adjacency_matrix.A, np.triu(adjacency_matrix.A, k=1)):
        raise ValueError(
            "The adjacency matrix should be upper triangular with zero diagonals."
        )
    return True


@optimod()
def maximum_weighted_independent_set(adjacency_matrix, weights, *, create_env):
    """Find a set of mutually non-adjacent vertices with maximum weighted sum.

    Parameters
    ----------
    adjacency_matrix : spmatrix
        The upper triangular adjacency matrix
    weights : ndarray
        Vertex weight array

    Returns
    -------
    Result
        A data class representing the maximum weighted independent set array
        and its weight
    """
    check_input(adjacency_matrix, weights)
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
        mwis = np.where(x.X >= 0.5)[0]
        return Result(mwis, sum(weights[mwis]))


def maximum_weighted_clique(adjacency_matrix, weights, **kwargs):
    """Find a set of fully connected vertices with maximum weighted sum.

    Parameters
    ----------
    adjacency_matrix : spmatrix
        The upper triangular adjacency matrix
    weights : ndarray
        Vertex weight array

    Returns
    -------
    Result
        A data class representing the maximum weighted independent set array
        and its weight
    """
    check_input(adjacency_matrix, weights)
    num_vertices, _ = adjacency_matrix.shape
    complement_matrix = (
        sp.triu(np.ones((num_vertices, num_vertices)), k=1) - adjacency_matrix
    )
    result = maximum_weighted_independent_set(complement_matrix, weights, **kwargs)
    return result
