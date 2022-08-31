import collections

import gurobipy as gp
from gurobipy import GRB
import scipy.sparse as sp


# Weighted graph as a sparse matrix.
row = [0, 1, 1, 2, 2, 3]
col = [3, 2, 3, 3, 4, 5]
data = [1, 1.2, 1.3, 1.4, 1, 1.2]
G = sp.coo_matrix((data, (row, col)), shape=(6, 6))

# Compute max matching.
m = gp.Model()

edges = list(zip(G.row, G.col))
# Assume G is bipartite, then this model is an LP
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

row, col, data = zip(*[(i, j, v.Obj) for (i, j), v in x.items() if v.X > 0.5])

matching = sp.coo_matrix((data, (row, col)), shape=G.shape)
