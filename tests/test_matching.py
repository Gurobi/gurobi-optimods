import unittest
from itertools import product
from random import Random

import numpy as np
import scipy.sparse as sp
from numpy.testing import assert_array_equal

from gurobi_optimods.matching import maximum_bipartite_matching, maximum_weighted_matching


def random_bipartite(n1, n2, m, seed):
    rng = Random(seed)
    edges = rng.sample(list(product(range(n1), range(n1, n1 + n2))), m)
    i, j = zip(*edges)
    data = np.ones(len(i))
    return sp.coo_matrix((data, (i, j)), shape=(n1 + n2, n1 + n2))


class TestBipartiteMatching(unittest.TestCase):
    def test_coo_matrix(self):
        G = random_bipartite(5, 5, 20, 0)
        matching = maximum_bipartite_matching(G)
        self.assertEqual(len(matching.data), 5)

    def test_csr_matrix(self):
        G = random_bipartite(5, 5, 20, 0).tocsr()
        matching = maximum_bipartite_matching(G)
        self.assertEqual(len(matching.data), 5)

    def test_not_bipartite(self):
        # Complete graph not bipartite, can recognise this by
        # fractional result
        data = [1, 1, 1]
        row = [0, 0, 1]
        col = [1, 2, 2]
        G = sp.coo_matrix((data, (row, col)), shape=(3, 3))
        with self.assertRaises(ValueError):
            maximum_bipartite_matching(G)


class TestWeightedMatching(unittest.TestCase):
    def test_coo_matrix(self):
        G = random_bipartite(5, 5, 20, 0)
        matching = maximum_weighted_matching(G)
        self.assertEqual(len(matching.data), 5)

    def test_csr_matrix(self):
        G = random_bipartite(5, 5, 20, 0).tocsr()
        matching = maximum_weighted_matching(G)
        self.assertEqual(len(matching.data), 5)

    def test_not_bipartite(self):
        # Complete graph not bipartite, can recognise this by
        # fractional result
        data = [1, 2, 3]
        row = [0, 0, 1]
        col = [1, 2, 2]
        G = sp.coo_matrix((data, (row, col)), shape=(3, 3))
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
