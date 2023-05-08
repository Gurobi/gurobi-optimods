import random
import unittest
from itertools import product

import numpy as np
import scipy.sparse as sp
from numpy.testing import assert_array_equal

from gurobi_optimods.general_matching import maximum_weighted_matching


def random_bipartite(n1, n2, p, seed):
    nodes1 = np.arange(n1)
    nodes2 = np.arange(n1, n1 + n2)
    rng = random.Random(seed)
    m = int(round(n1 * n2 * p))
    edges = rng.sample(list(product(nodes1, nodes2)), m)
    i, j = zip(*edges)
    data = np.ones(len(i))
    triu = sp.coo_array((data, (i, j)), shape=(n1 + n2, n1 + n2))
    return triu + triu.T, nodes1, nodes2


class TestGeneralMatching(unittest.TestCase):
    def test_coo_array(self):
        G, *_ = random_bipartite(5, 5, 0.8, 0)
        matching = maximum_weighted_matching(G)
        self.assertEqual(len(matching.data), 5)

    def test_csr_array(self):
        G, *_ = random_bipartite(5, 5, 0.8, 0)
        matching = maximum_weighted_matching(G.tocsr())
        self.assertEqual(len(matching.data), 5)

    def test_not_bipartite(self):
        # Complete graph not bipartite, general matching handles this
        data = [1, 2, 3]
        row = [0, 0, 1]
        col = [1, 2, 2]
        G = sp.coo_array((data, (row, col)), shape=(3, 3))
        matching = maximum_weighted_matching(G)
        self.assertEqual(len(matching.data), 1)
        expected = np.array(
            [
                [0, 0, 0],
                [0, 0, 3],
                [0, 0, 0],
            ]
        )
        assert_array_equal(matching.toarray(), expected)
