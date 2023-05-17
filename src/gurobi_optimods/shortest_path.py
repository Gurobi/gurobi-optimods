import logging
from typing import Optional, overload

import numpy as np
import pandas as pd
import scipy.sparse as sp

try:
    import networkx as nx
except ImportError:
    nx = None

from gurobi_optimods.network_util import solve_min_cost_flow
from gurobi_optimods.utils import optimod

logger = logging.getLogger(__name__)


@overload
def shortest_path(
    graph: sp.spmatrix,
    source: int,
    sink: int,
    silent: bool = False,
    logfile: Optional[str] = None,
) -> sp.spmatrix:
    ...


@overload
def shortest_path(
    graph: pd.DataFrame,
    source: int,
    sink: int,
    silent: bool = False,
    logfile: Optional[str] = None,
) -> pd.DataFrame:
    ...


if nx is not None:

    @overload
    def shortest_path(
        graph: nx.Graph,
        source: int,
        sink: int,
        silent: bool = False,
        logfile: Optional[str] = None,
    ) -> nx.Graph:
        ...


@optimod()
def shortest_path(graph, source: int, sink: int, *, create_env):
    """Solve the shortest path problem for a given graph.

    :param graph: A graph, specified either as a scipy.sparse adjacency matrix, networkx
        graph, or pandas dataframe
    :type graph: :class:`sp.sparray|nx.Graph|pd.DataFrame`
    :param source: The source node for the path.
    :type source: :class:`int`
    :param sink: The sink (or destination) node for the path.
    :type sink: :class:`int`
    :param silent: silent=True suppresses all console output (defaults to False)
    :type silent: bool
    :param logfile: Write all mod output to the given file path (defaults to None: no log)
    :type logfile: str
    :return: A subgraph of the original graph specifying the maximum matching
    :rtype: :class:`sp.sparray|nx.Graph|pd.DataFrame`
    """
    if isinstance(graph, sp.spmatrix):
        return _shortest_path_scipy(graph, source, sink, create_env)
    elif isinstance(graph, pd.DataFrame):
        return _shortest_path_pandas(graph, source, sink, create_env)
    elif nx is not None and isinstance(graph, nx.Graph):
        return _shortest_path_networkx(graph, source, sink, create_env)
    raise ValueError(f"Unknown graph type: {type(graph)}")


def _shortest_path_pandas(arc_data, source, sink, create_env):
    # Perform conversion from pd to ndarray
    # Arcs and attributes
    arcs = arc_data.index.to_numpy()
    from_arc = np.array([a[0] for a in arcs], dtype=np.int64)
    to_arc = np.array([a[1] for a in arcs], dtype=np.int64)

    cost = arc_data["cost"].to_numpy()

    capacity = np.ones(len(arcs), dtype=np.int64)

    # Nodes and attributes
    rep_nodes = np.concatenate((from_arc, to_arc), dtype=np.int64)
    # Get unique node labels in the order they appear in rep_nodes
    _, idx = np.unique(rep_nodes, return_index=True)
    nodes = rep_nodes[np.sort(idx)]

    demands = np.zeros(nodes.shape, dtype=object)
    demands[source] = -1.0
    demands[sink] = 1.0

    # Call solve_min_cost_flow using some data
    with create_env() as env:
        obj, flows = solve_min_cost_flow(env, from_arc, to_arc, capacity, cost, demands)

    return obj, pd.Series(flows, index=arc_data.index)


def _shortest_path_networkx(graph, source, sink, create_env):
    adjacency = nx.convert_matrix.to_scipy_sparse_array(graph)
    cost = nx.to_scipy_sparse_array(graph, weight="cost")
    adjacency.data = cost.data
    cost, sparray = _shortest_path_scipy(adjacency, source, sink, create_env)
    return cost, nx.convert_matrix.from_scipy_sparse_array(sparray)


def _shortest_path_scipy(G, source, sink, create_env):
    G = G.tocoo()

    from_arc = G.row
    to_arc = G.col

    costs = G.data
    capacities = np.ones(from_arc.shape, dtype=float)

    demands = np.zeros(G.shape[1], dtype=float)
    demands[source] = -1.0
    demands[sink] = 1.0

    with create_env() as env:
        cost, flows = solve_min_cost_flow(
            env, from_arc, to_arc, capacities, costs, demands
        )

    # Filter + create scipy output matrix
    select = flows > 0.5
    from_arc_result = from_arc[select]
    to_arc_result = to_arc[select]
    arg = (np.ones(from_arc_result.shape), (from_arc_result, to_arc_result))
    return cost, sp.coo_matrix(arg, dtype=int, shape=G.shape)
