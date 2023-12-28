"""
Maximum Weighted Independent Set/Clique
---------------------------------------
"""

from dataclasses import dataclass

import gurobipy as gp
import gurobipy_pandas as gppd
import numpy as np
import pandas as pd
import scipy.sparse as sp
from gurobipy import GRB

try:
    import networkx as nx
except ImportError:
    nx = None

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


def maximum_weighted_independent_set(graph, weights, **kwargs):
    """Find a set of mutually non-adjacent vertices with maximum weighted sum.

    Parameters
    ----------
    graph : spmatrix or Graph or DataFrame
        A graph, specified as a scipy.sparse upper triangular adjacency
        matrix with zero diagonals, networkx graph, or pandas dataframe.
        The pandas dataframe must include edge information in two columns
        named as ``"node1"`` and ``"node2"``.
    weights : ndarray or DataFrame
        Vertex weights. If ``graph`` is a scipy.sparse matrix or a
        networkx graph, ``weights`` must be a numpy array. If ``graph``
        is a pandas dataframe graph, ``weights`` must be a dataframe too.
        The weights dataframe must include the weight information in a column
        named ``"weights"`` and must be indexed by vertex number.

    Returns
    -------
    Result
        A data class representing the maximum weighted independent set array
        and its weight.
    """
    if sp.issparse(graph):
        return _maximum_weighted_independent_set_scipy(graph, weights, **kwargs)
    elif isinstance(graph, pd.DataFrame):
        return _maximum_weighted_independent_set_pandas(graph, weights, **kwargs)
    elif nx is not None and isinstance(graph, nx.Graph):
        return _maximum_weighted_independent_set_networkx(graph, weights, **kwargs)
    else:
        raise ValueError(f"Unknown graph type: {type(graph)}")


@optimod()
def _maximum_weighted_independent_set_scipy(adjacency_matrix, weights, *, create_env):
    """This implementation uses the gurobipy matrix friendly APIs which are well
    suited for the input data in scipy data structures."""
    with create_env() as env, gp.Model("mwis", env=env) as model:
        rows, cols = adjacency_matrix.tocoo().row, adjacency_matrix.tocoo().col
        num_vertices, num_edges = len(weights), len(rows)
        # x_i: 1 if vertex i is in the independent set and 0 otherwise
        x = model.addMVar(num_vertices, vtype=GRB.BINARY, name="x")
        # Maximize the sum of the vertex weights in the independent set
        model.setObjective(weights @ x, sense=GRB.MAXIMIZE)
        # Get the incident matrix from the adjacency matrix where
        # there is a column for each edge
        indices = []
        for i, j in zip(rows, cols):
            indices.extend([i, j])
        indptr = range(0, len(indices) + 2, 2)
        data = np.ones(2 * num_edges)
        A = sp.csc_array((data, indices, indptr), shape=(num_vertices, num_edges))
        # The independent set contains non-adjacent vertices
        model.addMConstr(
            A.T,
            x,
            GRB.LESS_EQUAL,
            np.ones(A.shape[1]),
            name="no_adjacent_vertices",
        )
        model.optimize()
        mwis = np.where(x.X >= 0.5)[0]
        return Result(mwis, sum(weights[mwis]))


@optimod()
def _maximum_weighted_independent_set_pandas(frame, weights, *, create_env):
    """This implementation uses the gurobipy-pandas APIs which are well
    suited for the input data in pandas dataframes structures."""
    with create_env() as env, gp.Model("mwis", env=env) as model:
        # x_i: 1 if vertex i is in the independent set and 0 otherwise
        x = gppd.add_vars(model, weights, name="x", vtype=GRB.BINARY)
        # Maximize the sum of the vertex weights in the independent set
        model.setObjective((x * weights["weights"]).sum(), sense=GRB.MAXIMIZE)

        # The independent set contains non-adjacent vertices
        def f(node1, node2):
            return x[node1] + x[node2]

        if len(frame) > 0:
            gppd.add_constrs(
                model,
                frame.apply(lambda x: f(x["node1"], x["node2"]), axis=1),
                GRB.LESS_EQUAL,
                1,
                name="no_adjacent_vertices",
            )
        model.optimize()
        (mwis,) = np.where(x.gppd.X >= 0.5)
        return Result(mwis, weights["weights"].iloc[mwis].sum())


@optimod()
def _maximum_weighted_independent_set_networkx(graph, weights, *, create_env):
    """This implementation uses the gurobipy term-based APIs which are well
    suited for the input data in networkx data structures."""
    with create_env() as env, gp.Model("mwis", env=env) as model:
        num_nodes, edges = len(weights), graph.edges
        # x_i: 1 if vertex i is in the independent set and 0 otherwise
        x = model.addVars(num_nodes, vtype=GRB.BINARY, name="x")
        # Maximize the sum of the vertex weights in the independent set
        model.setObjective(
            gp.quicksum(x[node] * weights[node] for node in range(num_nodes)),
            sense=GRB.MAXIMIZE,
        )
        # The independent set contains non-adjacent vertices
        model.addConstrs(
            (x[node1] + x[node2] <= 1 for (node1, node2) in edges),
            name="no_adjacent_vertices",
        )
        model.optimize()
        (mwis,) = np.where(model.getAttr("X", model.getVars()))
        return Result(mwis, sum(weights[mwis]))


def maximum_weighted_clique(graph, weights, **kwargs):
    """Find a set of fully connected vertices with maximum weighted sum.

    Parameters
    ----------
    graph : spmatrix or Graph or DataFrame
        A graph, specified as a scipy.sparse upper triangular adjacency
        matrix with zero diagonals, networkx graph, or pandas dataframe.
        The pandas dataframe must include edge information in two columns
        named as ``"node1"`` and ``"node2"``.
    weights : ndarray of DataFrame
        Vertex weights. If ``graph`` is a scipy.sparse matrix or a
        networkx graph, ``weightss` must be a numpy array. If ``graph``
        is a pandas dataframe graph, ``weights`` must be a dataframe too.
        The weights dataframe must include the weight information in a column
        named ``"weights"`` and must be indexed by vertex number.

    Returns
    -------
    Result
        A data class representing the maximum weighted clique array
        and its weight.
    """
    if sp.issparse(graph):
        num_vertices, _ = graph.shape
        complement_matrix = sp.triu(np.ones((num_vertices, num_vertices)), k=1) - graph
        return _maximum_weighted_independent_set_scipy(
            complement_matrix, weights, **kwargs
        )
    elif isinstance(graph, pd.DataFrame):
        num_vertices = len(weights)
        data = (
            [node1, node2]
            for node1 in range(num_vertices)
            for node2 in range(node1 + 1, num_vertices)
            if (node1, node2) not in zip(graph["node1"], graph["node2"])
        )
        complement_frame = pd.DataFrame(data, columns=["node1", "node2"])
        return _maximum_weighted_independent_set_pandas(
            complement_frame, weights, **kwargs
        )
    elif nx is not None and isinstance(graph, nx.Graph):
        return _maximum_weighted_independent_set_networkx(
            nx.complement(graph), weights, **kwargs
        )
    else:
        raise ValueError(f"Unknown graph type: {type(graph)}")
