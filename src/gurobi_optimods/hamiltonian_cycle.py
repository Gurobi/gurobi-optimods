"""
Bottleneck Hamiltonian Cycle
----------------------------
"""

import logging
from dataclasses import dataclass
from typing import Any

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


@dataclass
class HamiltonianCycleResult:
    """
    Data class representing a directed Hamiltonian cycle and the largest
    arc weight it uses.

    Attributes
    ----------
    cycle : list
        Nodes of the graph in visiting order. Every node appears exactly once,
        and the cycle is closed by returning from the last node to the first.
    weights : list
        Weight of each arc of the cycle, aligned with ``cycle`` so that
        ``weights[k]`` is the weight of the arc from ``cycle[k]`` to
        ``cycle[k + 1]``, wrapping around at the end.
    bottleneck : Any
        The largest value in ``weights``. This is the quantity minimized.
    """

    cycle: list
    weights: list
    bottleneck: Any


@optimod()
def min_bottleneck_hamiltonian_cycle(graph, *, create_env):
    """Find a directed Hamiltonian cycle minimizing the largest arc weight.

    A Hamiltonian cycle visits every node of the graph exactly once before
    returning to its starting point. Among all such cycles, this Mod returns
    one whose largest arc weight is as small as possible. This is known as the
    bottleneck (or min-max) Hamiltonian cycle problem.

    Parameters
    ----------
    graph : spmatrix or sparray or Graph or DiGraph or DataFrame
        A directed graph, specified as a scipy.sparse adjacency matrix whose
        stored entries are the arc weights, a networkx graph whose edges carry
        a ``"weight"`` attribute, or a pandas dataframe indexed by
        ``("source", "target")`` with a ``"weight"`` column. Weights are
        optional everywhere and default to one, in which case the Mod returns
        any Hamiltonian cycle. An undirected networkx graph is read as a
        directed graph in which every edge may be traversed either way.

    Returns
    -------
    HamiltonianCycleResult
        A data class containing the cycle, the weight of each of its arcs, and
        the largest of those weights.

    Raises
    ------
    ValueError
        If the graph has fewer than three nodes, if some node has no incoming
        or no outgoing arc, or if the graph contains no Hamiltonian cycle.
    """
    if sp.issparse(graph):
        nodes, tails, heads, weights = _extract_scipy(graph)
    elif isinstance(graph, pd.DataFrame):
        nodes, tails, heads, weights = _extract_pandas(graph)
    elif nx is not None and isinstance(graph, nx.Graph):
        nodes, tails, heads, weights = _extract_networkx(graph)
    else:
        raise ValueError(f"Unknown graph type: {type(graph)}")

    return _solve_bottleneck(nodes, tails, heads, weights, create_env)


def _extract_scipy(graph):
    """Stored entries of the sparse matrix define the arcs, and their values
    are the arc weights."""
    matrix = graph.tocoo()
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"Adjacency matrix must be square, got {matrix.shape}")
    nodes = list(range(matrix.shape[0]))
    return nodes, matrix.row.astype(int), matrix.col.astype(int), matrix.data


def _extract_pandas(frame):
    """Arcs come from the ("source", "target") index, or from columns of those
    names, and weights from an optional "weight" column."""
    if {"source", "target"}.issubset(frame.columns):
        source = frame["source"].to_numpy()
        target = frame["target"].to_numpy()
    elif isinstance(frame.index, pd.MultiIndex) and frame.index.nlevels == 2:
        source = frame.index.get_level_values(0).to_numpy()
        target = frame.index.get_level_values(1).to_numpy()
    else:
        raise ValueError(
            "Dataframe must be indexed by ('source', 'target'), or have "
            "'source' and 'target' columns"
        )

    if "weight" in frame.columns:
        weights = frame["weight"].to_numpy()
    else:
        weights = np.ones(len(frame), dtype=int)

    nodes = list(pd.unique(np.concatenate([source, target])))
    index = {node: i for i, node in enumerate(nodes)}
    tails = np.array([index[node] for node in source], dtype=int)
    heads = np.array([index[node] for node in target], dtype=int)
    return nodes, tails, heads, weights


def _extract_networkx(graph):
    """Edges carry their weight in a "weight" attribute. Edges of an
    undirected graph may be traversed in either direction."""
    nodes = list(graph.nodes)
    index = {node: i for i, node in enumerate(nodes)}
    directed = graph.is_directed()

    tails, heads, weights = [], [], []
    for source, target, data in graph.edges(data=True):
        weight = data.get("weight", 1)
        tails.append(index[source])
        heads.append(index[target])
        weights.append(weight)
        if not directed:
            tails.append(index[target])
            heads.append(index[source])
            weights.append(weight)

    return (
        nodes,
        np.array(tails, dtype=int),
        np.array(heads, dtype=int),
        np.array(weights),
    )


def _to_python(value):
    """Return numpy scalars in their most natural python form, so that results
    read the same way as the data the user handed in."""
    if isinstance(value, np.datetime64):
        return pd.Timestamp(value)
    if isinstance(value, np.generic):
        return value.item()
    return value


def _prepare(nodes, tails, heads, weights):
    """Clean up the arc list and reduce the weights to ranks.

    Only the ordering of the weights matters to a bottleneck objective, so the
    model is built on the ranks of the distinct weight values. This keeps the
    objective integral and small, and means datetimes need no special casing.
    """
    num_nodes = len(nodes)
    if num_nodes < 3:
        raise ValueError(
            f"A directed Hamiltonian cycle needs at least 3 nodes, got {num_nodes}"
        )

    tails = np.asarray(tails, dtype=int)
    heads = np.asarray(heads, dtype=int)
    weights = np.asarray(weights)

    # Self loops can never appear in a Hamiltonian cycle
    keep = tails != heads
    tails, heads, weights = tails[keep], heads[keep], weights[keep]
    if len(tails) == 0:
        raise ValueError("Graph has no arcs, so it contains no Hamiltonian cycle")

    # Replace weights by their rank among the distinct values present
    values, ranks = np.unique(weights, return_inverse=True)
    ranks = ranks.reshape(-1)

    # Collapse parallel arcs, keeping the smallest weight for each ordered
    # pair: one arc is enough, and the earliest is never worse.
    order = np.lexsort((ranks, heads, tails))
    tails, heads, ranks = tails[order], heads[order], ranks[order]
    first = np.ones(len(tails), dtype=bool)
    first[1:] = (tails[1:] != tails[:-1]) | (heads[1:] != heads[:-1])
    tails, heads, ranks = tails[first], heads[first], ranks[first]

    # Every node needs one arc in and one arc out, so a node missing either
    # rules out a Hamiltonian cycle. Say so, rather than reporting infeasible.
    for label, direction in ((tails, "outgoing"), (heads, "incoming")):
        counts = np.bincount(label, minlength=num_nodes)
        missing = np.flatnonzero(counts == 0)
        if missing.size:
            raise ValueError(
                f"Node {_to_python(nodes[missing[0]])!r} has no {direction} arcs, "
                f"so the graph contains no Hamiltonian cycle"
            )

    # Each node must use one of its outgoing arcs and one of its incoming
    # arcs, so the bottleneck is at least the cheapest of either, for every
    # node. This bound is free and is often strong.
    large = np.iinfo(ranks.dtype).max
    out_min = np.full(num_nodes, large, dtype=ranks.dtype)
    in_min = np.full(num_nodes, large, dtype=ranks.dtype)
    np.minimum.at(out_min, tails, ranks)
    np.minimum.at(in_min, heads, ranks)
    lower_bound = int(max(out_min.max(), in_min.max()))

    return tails, heads, ranks, values, lower_bound


def _subtour_callback(model, where):
    """Reject solutions made up of several disjoint cycles.

    The degree constraints already force every node to have one arc in and one
    arc out, so an integer solution is always a disjoint union of cycles. If
    there is more than one, each is cut off.
    """
    if where != GRB.Callback.MIPSOL:
        return

    values = model.cbGetSolution(model._x)
    successor = {i: j for (i, j) in model._arcs if values[i, j] > 0.5}

    unvisited = set(successor)
    while unvisited:
        component = []
        node = next(iter(unvisited))
        while node in unvisited:
            unvisited.remove(node)
            component.append(node)
            node = successor[node]

        if len(component) < model._num_nodes:
            inside = set(component)
            model.cbLazy(
                gp.quicksum(
                    model._x[i, j]
                    for (i, j) in model._arcs
                    if i in inside and j in inside
                )
                <= len(component) - 1
            )


def _solve_bottleneck(nodes, tails, heads, weights, create_env):
    tails, heads, ranks, values, lower_bound = _prepare(nodes, tails, heads, weights)
    num_nodes, num_arcs = len(nodes), len(tails)
    arcs = list(zip(tails.tolist(), heads.tolist()))

    logger.info(
        f"Solving bottleneck Hamiltonian cycle over {num_nodes} nodes "
        f"and {num_arcs} arcs"
    )

    with create_env() as env, gp.Model("hamiltonian_cycle", env=env) as model:
        # x[i, j]: 1 if the cycle travels directly from node i to node j
        x = model.addVars(arcs, vtype=GRB.BINARY, name="x")
        # z: the largest weight used by the cycle, which we minimize
        z = model.addVar(
            vtype=GRB.INTEGER, lb=lower_bound, ub=len(values) - 1, name="z"
        )
        model.setObjective(z, GRB.MINIMIZE)

        # Every node is entered exactly once and left exactly once
        model.addConstrs((x.sum(i, "*") == 1 for i in range(num_nodes)), name="out")
        model.addConstrs((x.sum("*", j) == 1 for j in range(num_nodes)), name="in")

        # z bounds every arc the cycle uses. Arcs at or below the lower bound
        # cannot raise z, so they need no constraint.
        for (i, j), rank in zip(arcs, ranks.tolist()):
            if rank > lower_bound:
                model.addConstr(z >= rank * x[i, j], name=f"bottleneck[{i},{j}]")

        # Subtour elimination constraints are added lazily as they are needed
        model._x = x
        model._arcs = arcs
        model._num_nodes = num_nodes
        model.Params.LazyConstraints = 1
        model.optimize(_subtour_callback)

        if model.Status in (GRB.INFEASIBLE, GRB.INF_OR_UNBD):
            raise ValueError("Graph contains no Hamiltonian cycle")
        if model.SolCount == 0:
            raise ValueError(
                "No Hamiltonian cycle was found before the solve terminated"
            )
        if model.Status != GRB.OPTIMAL:
            logger.info(
                "Solve terminated early; returning the best cycle found, which "
                "may not minimize the bottleneck"
            )

        solution = model.getAttr("X", x)

    # Walk the successor function to read the cycle off in visiting order
    successor = {i: j for (i, j) in arcs if solution[i, j] > 0.5}
    order = [0]
    while successor[order[-1]] != 0:
        order.append(successor[order[-1]])

    rank_of = dict(zip(arcs, ranks.tolist()))
    cycle_weights = [
        _to_python(values[rank_of[(order[k], order[(k + 1) % num_nodes])]])
        for k in range(num_nodes)
    ]
    bottleneck = max(cycle_weights)
    logger.info(f"Bottleneck weight of the cycle found: {bottleneck}")

    return HamiltonianCycleResult(
        cycle=[_to_python(nodes[i]) for i in order],
        weights=cycle_weights,
        bottleneck=bottleneck,
    )
