import logging
from typing import Optional, overload

import gurobipy as gp
import numpy as np
import pandas as pd
import scipy.sparse as sp
from gurobipy import GRB

try:
    import networkx as nx
except ImportError:
    nx = None

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


def solve_min_cost_flow(
    env: gp.Env,
    edge_source: np.ndarray,
    edge_target: np.ndarray,
    capacity: np.ndarray,
    cost: np.ndarray,
    demand: np.ndarray,
):
    """Solve a min-cost flow model for the given network.

    Formulated as min-cost flow (A, b, c, u) using the matrix-friendly API, where
    A = incidence matrix, b = demand, c = cost, u = capacity::

        min   c^T x
        s.t.  A x = b
              0 <= x <= u

    :param env: Gurobi environment to use for model
    :type env: :class:`gurobipy.Env`
    :param edge_source: Index array; ``edge_source[i]`` is the source node
        :math:`j` of edge :math:`i`
    :type edge_source: :class:`numpy.ndarray`
    :param edge_target: Index array; ``edge_target[i]`` is the source node
        :math:`j` of edge :math:`i`
    :type edge_target: :class:`numpy.ndarray`
    :param capacity: Capacity array; ``capacity[i]`` is the capacity of
        edge :math:`i`
    :type capacity: :class:`numpy.ndarray`
    :param cost: Cost array; ``cost[i]`` is the cost of edge :math:`i` (to
        be minimized)
    :type cost: :class:`numpy.ndarray`
    :param demand: Demand array; ``demand[j]`` is the net demand of node
        :math:`j`
    :type demand: :class:`numpy.ndarray`
    """

    # Create incidence matrix from edge lists.
    indices = np.column_stack((edge_source, edge_target)).reshape(-1, order="C")
    indptr = np.arange(0, 2 * edge_source.shape[0] + 2, 2)
    ones = np.ones(edge_source.shape)
    data = np.column_stack((ones * -1.0, ones)).reshape(-1, order="C")
    A = sp.csc_array((data, indices, indptr))

    logger.info("Solving min-cost flow with {0} nodes and {1} edges".format(*A.shape))

    # Solve model with gurobi, return cost and flows
    with gp.Model(env=env) as m:
        m.ModelSense = GRB.MINIMIZE
        x = m.addMVar(A.shape[1], lb=0, ub=capacity, obj=cost, name="x")
        m.addMConstr(A, x, GRB.EQUAL, demand)
        m.optimize()
        if m.Status == GRB.INFEASIBLE:
            raise ValueError("Unsatisfiable flows")
        return m.ObjVal, x.X
