import collections

import numpy as np
import gurobipy as gp
import scipy.sparse as sp


def maximum_matching(G):
    """Return a subgraph which is the maximum matching of G."""

    m = gp.Model()

    G = G.tocoo()
    edges = list(zip(G.row, G.col))
    x = m.addVars(edges, name="x", vtype=gp.GRB.BINARY)

    clashes = collections.defaultdict(set)
    for edge in edges:
        clashes[edge[0]].add(edge)
        clashes[edge[1]].add(edge)

    for edges in clashes.values():
        m.addConstr(gp.quicksum(x[edge] for edge in edges) <= 1)

    m.setObjective(x.sum(), sense=gp.GRB.MAXIMIZE)
    m.optimize()

    selected = {edge for edge, v in x.items() if v.X > 0.5}
    i, j = zip(*selected)
    data = np.ones(len(i))

    return sp.coo_matrix((data, (i, j)), shape=G.shape)
