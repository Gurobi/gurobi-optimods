import collections

import numpy as np
import gurobipy as gp
from gurobipy import GRB
import scipy.sparse as sp


def maximum_bipartite_matching(G):
    """Return a subgraph which is the maximum cardinality matching of
    the bipartite graph G.

    :param G: Adjacency matrix of an unweighted bipartite graph.
    :type G: :class:`sp.spmatrix`
    :return: Adjacency matrix of the maximum bipartite matching subgraph
    :rtype: :class:`sp.spmatrix`
    """

    with gp.Env() as env, gp.Model(env=env) as m:

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

        return sp.coo_matrix((data, (i, j)), shape=G.shape)


def maximum_weighted_matching(G):
    """Return a subgraph which is the maximum weighted matching of G.

    :param G: Adjacency matrix of a unweighted graph.
    :type G: :class:`sp.spmatrix`
    :return: Adjacency matrix of the maximum weighted matching subgraph
    :rtype: :class:`sp.spmatrix`
    """

    with gp.Env() as env, gp.Model(env=env) as m:

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

        return sp.coo_matrix((data, (row, col)), shape=G.shape)
