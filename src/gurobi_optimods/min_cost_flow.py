import collections
import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd
import gurobipy as gp
from gurobipy import GRB
import scipy.sparse as sp

try:
    import networkx as nx
except ImportError:
    nx = None

from gurobi_optimods.network_util import solve_min_cost_flow
from gurobi_optimods.utils import optimod

logger = logging.getLogger(__name__)


@optimod()
def min_cost_flow(arc_data: pd.DataFrame, demand_data: pd.DataFrame, *, create_env):
    """Solve the minimum cost flow problem for a given graph.

    The inputs adhere to the following structure:

    .. code-block:: python

        arc_data = pd.DataFrame(
            [
                {"from": 0, "to": 1, "capacity": 16, "cost": 0},
                {"from": 1, "to": 2, "capacity": 10, "cost": 0},
            ]
        ).set_index(["from", "to"])

        demand_data = pd.DataFrame(
            [{"node": 0, "demand": -1}, {"node": 2, "demand": 1}]
        ).set_index("node")

    :param arc_data: DataFrame with graph and respective attributes. These must include ``"from"``, ``"to"`` nodes used as index, ``"capacity"``, and ``"cost"``.
    :type arc_data: :class:`pd.DataFrame`
    :param demand_data: DataFrame with node demand information. These must include indexed by `"node"`, and include the `"demand"`. This value can be positive (requesting flow) or negative (supplying flow).
    :type demand_data: :class:`pd.DataFrame`
    :param silent: Optional. Boolean with whether output should be printed.
    :type silent: :class:`bool`
    :param logfile: Optional. String with file path with logger and Gurobi.
    :type logfile: :class:`str`
    :return: Cost of the minimum cost flow.
    :rtype: :class:`float`
    :return: DataFrame with the flow for each edge.
    :rtype: :class:`pd.Series`
    """
    # Perform conversion from pd to ndarray
    # Arcs and attributes
    arcs = arc_data.index.to_numpy()
    from_arc = np.array([a[0] for a in arcs], dtype="object")
    to_arc = np.array([a[1] for a in arcs], dtype="object")
    capacity = arc_data["capacity"].to_numpy()
    cost = arc_data["cost"].to_numpy()

    # Nodes and attributes
    rep_nodes = np.concatenate((from_arc, to_arc), dtype="object")
    # Get unique node labels in the order they appear in rep_nodes
    _, idx = np.unique(rep_nodes, return_index=True)
    nodes = rep_nodes[np.sort(idx)]

    demand_nodes = demand_data.index.to_numpy()
    _demands = demand_data["demand"].to_numpy()
    demands = np.zeros(nodes.shape, dtype=object)
    # Assign demand in the same order as the nodes array
    demand_nodes_indeces = np.where(np.isin(nodes, demand_nodes, assume_unique=True))[0]
    for i, n in enumerate(demand_nodes_indeces):
        demands[n] = _demands[i]

    # Call solve_min_cost_flow using some data
    with create_env() as env:
        obj, flows = solve_min_cost_flow(env, from_arc, to_arc, capacity, cost, demands)

    return obj, pd.Series(flows, index=arc_data.index)


@optimod()
def min_cost_flow_scipy(
    G: sp.spmatrix,
    capacities: sp.spmatrix,
    costs: sp.spmatrix,
    demands: np.ndarray,
    *,
    create_env,
):
    """Solve the minimum cost flow problem for a given graph.

    :param G: Adjacency matrix of the graph.
    :type G: :class:`sp.sparray`
    :param capacities: Matrix containing capacities for each edge.
    :type capacities: :class:`sp.sparray`
    :param costs: Matrix containing costs for each edge.
    :type costs: :class:`sp.sparray`
    :param demands: Array containing the demand for each node.
    :type demands: :class:`np.ndarray`
    :param silent: Optional. Boolean with whether output should be printed.
    :type silent: :class:`bool`
    :param logfile: Optional. String with file path with logger and Gurobi.
    :type logfile: :class:`str`
    :return: Cost of the minimum cost flow.
    :rtype: :class:`float`
    :return: Adjacency matrix with flow in the solution
    :rtype: :class:`sp.sparray`
    """
    G = G.tocoo()

    from_arc = G.row
    to_arc = G.col

    capacities = capacities.tocoo()
    capacities = capacities.data

    costs = costs.tocoo()
    costs = costs.data

    with create_env() as env:
        cost, flows = solve_min_cost_flow(
            env, from_arc, to_arc, capacities, costs, demands
        )
    # Filter + create scipy output matrix
    select = flows > 0.5
    from_arc_result = from_arc[select]
    to_arc_result = to_arc[select]
    arg = (flows[select], (from_arc_result, to_arc_result))
    return cost, sp.coo_matrix(arg, dtype=float, shape=G.shape)


@optimod()
def min_cost_flow_networkx(G, *, create_env):
    """Solve the minimum cost flow problem for a given graph.

    Note: We assume the networkx input graph node labels are all integers.

    :param G: Graph with edge attributes ``capacity`` and ``cost``, as well as node attributes ``demand``.
    :type G: :class:`nx.DiGraph`
    :param silent: Optional. Boolean with whether output should be printed.
    :type silent: :class:`bool`
    :param logfile: Optional. String with file path with logger and Gurobi.
    :type logfile: :class:`str`
    :return: Cost of the minimum cost flow.
    :rtype: :class:`float`
    :return: Dictionary indexed by edges with non-zero flow in the solution.
    :rtype: :class:`dict`
    """
    # Perform conversion from nx to ndarray
    edges = list(G.edges(data=True))  # ensure ordering

    from_arc = np.array([e[0] for e in edges])
    to_arc = np.array([e[1] for e in edges])
    capacity = np.array([e[2]["capacity"] for e in edges])
    cost = np.array([e[2]["cost"] for e in edges])

    nodes = list(G.nodes(data=True))  # ensure ordering
    demands = np.zeros(len(nodes))
    for i, v in enumerate(nodes):
        n, d = v  # unpack tuple
        if d:
            demands[i] = d["demand"]

    with create_env() as env:
        obj, flows = solve_min_cost_flow(env, from_arc, to_arc, capacity, cost, demands)

    # Convert numpy sol back to networkx
    select = flows > 0.5
    from_arc_result = from_arc[select]
    to_arc_result = to_arc[select]
    return obj, {
        (f, t): v for f, t, v in zip(from_arc_result, to_arc_result, flows[select])
    }
