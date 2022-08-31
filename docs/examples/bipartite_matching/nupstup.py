import scipy.sparse as sp

from nupstup.matching import maximum_bipartite_matching

# Bipartite graph as a sparse matrix.
row = [0, 3, 4, 0, 1, 3]
col = [7, 5, 5, 6, 6, 7]
data = [1, 1, 1, 1, 1, 1]
G = sp.coo_matrix((data, (row, col)))

# Compute max matching.
matching = maximum_bipartite_matching(G)
