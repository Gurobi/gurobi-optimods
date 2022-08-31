import scipy.sparse as sp

from nupstup.matching import maximum_weighted_matching

# Weighted graph as a sparse matrix.
row = [0, 1, 1, 2, 2, 3]
col = [3, 2, 3, 3, 4, 5]
data = [1, 1.2, 1.3, 1.4, 1, 1.2]
G = sp.coo_matrix((data, (row, col)), shape=(6, 6))

# Compute max matching.
matching = maximum_weighted_matching(G)
