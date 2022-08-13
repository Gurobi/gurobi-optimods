import unittest
from itertools import product
from random import Random

import numpy as np
import scipy.sparse as sp

from nupstup.matching import maximum_matching


def random_bipartite(n1, n2, m, seed):
    rng = Random(seed)
    edges = rng.sample(list(product(range(n1), range(n1, n1 + n2))), m)
    i, j = zip(*edges)
    data = np.ones(len(i))
    return sp.coo_matrix((data, (i, j)), shape=(n1 + n2, n1 + n2))


class TestMatching(unittest.TestCase):
    def test_coo_matrix(self):
        G = random_bipartite(5, 5, 20, 0)
        matching = maximum_matching(G)
        self.assertEqual(len(matching.data), 5)

    def test_csr_matrix(self):
        G = random_bipartite(5, 5, 20, 0).tocsr()
        matching = maximum_matching(G)
        self.assertEqual(len(matching.data), 5)
