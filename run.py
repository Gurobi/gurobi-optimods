import networkx as nx
import numpy as np
import scipy.sparse as sp

from gurobi_optimods.matching import maximum_bipartite_matching


n1, n2, p, seed = 3, 5, 0.5, 2394
nodes1 = np.arange(n1)
nodes2 = np.arange(n1, n1 + n2)
graph = nx.bipartite.random_graph(n1, n2, p, seed, directed=False)
adjacency = nx.to_scipy_sparse_array(graph)
print("=== Input graph ===")
print(sp.triu(adjacency).tocoo())

print("=== Loud mod ===")
matching = maximum_bipartite_matching(
    adjacency, nodes1, nodes2, silent=False, logfile="matching.log"
)
print(sp.triu(matching).tocoo())

print("=== Quiet mod ===")
matching = maximum_bipartite_matching(adjacency, nodes1, nodes2, silent=True)
print(sp.triu(matching).tocoo())
