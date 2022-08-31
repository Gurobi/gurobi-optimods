import collections

import gurobipy as gp
import scipy.sparse as sp


# Bipartite graph as a sparse matrix.
row = [0, 3, 4, 0, 1, 3]
col = [7, 5, 5, 6, 6, 7]
data = [1, 1, 1, 1, 1, 1]
G = sp.coo_matrix((data, (row, col)))

# Compute max matching.
m = gp.Model()

edges = list(zip(G.row, G.col))
# Assume G is bipartite, then this model is an LP
x = m.addVars(edges, name="x")

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
data = [1 for _ in i]

matching = sp.coo_matrix((data, (i, j)), shape=G.shape)
