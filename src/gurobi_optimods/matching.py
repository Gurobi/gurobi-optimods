import collections
import logging

import numpy as np
import gurobipy as gp
from gurobipy import GRB
import scipy.sparse as sp

from gurobi_optimods.utils import optimod

logger = logging.getLogger(__name__)


@optimod(mod_logger=logger)
def maximum_bipartite_matching(G, *, create_env):
    """Return a subgraph which is the maximum cardinality matching of
    the bipartite graph G.

    :param G: Adjacency matrix of an unweighted bipartite graph.
    :type G: :class:`sp.sparray`
    :return: Adjacency matrix of the maximum bipartite matching subgraph
    :rtype: :class:`sp.sparray`
    """

    logger.info(
        f"Solving bipartite matching model with {G.shape[0]} nodes and "
        f"{int(G.nnz/2)} edges"
    )

    with create_env() as env, gp.Model(env=env) as m:
        G = G.tocoo()
        edges = list(zip(G.row, G.col))
        # Continuous ok for bipartite case
        x = m.addVars(edges, name="x")

        clashes = collections.defaultdict(set)
        for edge in edges:
            clashes[edge[0]].add(edge)
            clashes[edge[1]].add(edge)

        for edgepair in clashes.values():
            m.addConstr(gp.quicksum(x[edge] for edge in edgepair) <= 1)

        m.setObjective(x.sum(), sense=GRB.MAXIMIZE)
        m.optimize()

        # FIXME appropriate tolerance?
        if any(abs(v.X - round(v.X)) > 1e-4 for v in x.values()):
            raise ValueError("Input graph not bipartite")

        selected = {edge for edge, v in x.items() if v.X > 0.5}
        i, j = zip(*selected)
        data = np.ones(len(i))

        logger.info(f"Max bipartite matching has {data.shape[0]} edges")

        return sp.coo_array((data, (i, j)), shape=G.shape)


def solve_min_cost_flow(env, from_arc, to_arc, capacity, cost):
    # Formulate as min-cost flow (A, b, c, u)
    #
    #   min  c^T x
    #  s.t.  A x = b
    #        0 <= x <= u
    #
    s = from_arc.shape[0]
    i = np.concatenate([from_arc, to_arc])
    j = np.concatenate([np.arange(s), np.arange(s)])
    data = np.ones(s * 2, dtype=float)
    data[s:] *= -1.0
    A = sp.coo_matrix((data, (i, j))).tocsr()
    b = np.zeros(A.shape[0])
    u = capacity
    c = cost
    names = np.array([f"x[{i},{j}]" for i, j in zip(from_arc, to_arc)])

    # Solve model with gurobi, return flows
    with gp.Model(env=env) as m:
        m.ModelSense = GRB.MINIMIZE
        x = m.addMVar(s, lb=0, ub=u, obj=c, name=names)
        m.addMConstr(A, x, GRB.EQUAL, b)
        m.optimize()
        return x.X


def maximum_bipartite_matching_flow(adjacency, nodes1, nodes2):
    G_nodes = adjacency.shape[0]
    source, sink = G_nodes, G_nodes + 1
    N_nodes = G_nodes + 2
    G = sp.triu(adjacency.tocoo())

    logger.info(
        f"Solving maximum matching n1={nodes1.shape[0]} "
        f"n2={nodes2.shape[0]} |E|={G.data.shape[0]}"
    )

    # Build network:
    #   source -> nodes1 (complete)
    #   nodes1 -> nodes2 (adjacency)
    #   nodes2 -> sink (complete)
    #   sink -> source
    from_arc = np.concatenate([[source] * nodes1.shape[0], G.row, nodes2, [sink]])
    to_arc = np.concatenate([nodes1, G.col, [sink] * nodes2.shape[0], [source]])
    capacity = np.ones(from_arc.shape, dtype=float)
    capacity[-1] = GRB.INFINITY
    cost = np.zeros(from_arc.shape, dtype=float)
    cost[-1] = -1.0
    logger.info(
        "Maximum matching formulated as min-cost flow "
        f"with {N_nodes} nodes and {from_arc.shape[0]} arcs"
    )

    # Solve min-cost flow problem
    with gp.Env() as env:
        flows = solve_min_cost_flow(env, from_arc, to_arc, capacity, cost)

    # Choose the arcs corresponding to the original graph with non-zero
    # flow. Note that the last var is the sink->source connection (drop it).
    select = (flows > 0.5) & (from_arc != source) & (to_arc != sink)
    from_arc_result = from_arc[select][:-1]
    to_arc_result = to_arc[select][:-1]

    logger.info(f"Done: max bipartite matching has {from_arc_result.shape[0]} edges")

    # Return undirected, unweighted adjacency matrix
    arg = (np.ones(from_arc_result.shape), (from_arc_result, to_arc_result))
    matching = sp.coo_matrix(arg, dtype=float, shape=G.shape)
    return matching + matching.T


@optimod(mod_logger=logger)
def maximum_weighted_matching(G, *, create_env):
    """Return a subgraph which is the maximum weighted matching of G.

    :param G: Adjacency matrix of a unweighted graph.
    :type G: :class:`sp.sparray`
    :return: Adjacency matrix of the maximum weighted matching subgraph
    :rtype: :class:`sp.sparray`
    """

    logger.info(
        f"Solving weighted matching model with {G.shape[0]} nodes and "
        f"{int(G.nnz/2)} edges"
    )

    with create_env() as env, gp.Model(env=env) as m:
        G = G.tocoo()
        edges = list(zip(G.row, G.col))
        x = m.addVars(edges, name="x", vtype=GRB.BINARY)

        clashes = collections.defaultdict(set)
        for edge in edges:
            clashes[edge[0]].add(edge)
            clashes[edge[1]].add(edge)

        for edgepair in clashes.values():
            m.addConstr(gp.quicksum(x[edge] for edge in edgepair) <= 1)

        weights = dict(zip(edges, G.data))
        m.setObjective(x.prod(weights), sense=GRB.MAXIMIZE)
        m.optimize()

        # FIXME appropriate tolerance?
        if any(abs(v.X - round(v.X)) > 1e-4 for v in x.values()):
            raise ValueError("Input graph not bipartite")

        row, col, data = zip(*[(i, j, v.Obj) for (i, j), v in x.items() if v.X > 0.5])

        logger.info(f"Max weighted matching has {len(data)} edges")

        return sp.coo_array((data, (row, col)), shape=G.shape)
