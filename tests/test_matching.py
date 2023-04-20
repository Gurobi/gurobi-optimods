import unittest

import networkx as nx
import numpy as np
import scipy.sparse as sp
from numpy.testing import assert_array_equal

from gurobi_optimods.matching import (
    maximum_bipartite_matching,
    maximum_weighted_matching,
)


def random_bipartite(n1, n2, p, seed):
    nodes1 = np.arange(n1)
    nodes2 = np.arange(n1, n1 + n2)
    graph = nx.bipartite.random_graph(n1, n2, p, seed, directed=False)
    adjacency = nx.to_scipy_sparse_array(graph)
    return adjacency, nodes1, nodes2


class TestBipartiteMatching(unittest.TestCase):
    def test_coo_array(self):
        G, *_ = random_bipartite(5, 5, 0.5, 0)
        matching = maximum_bipartite_matching(G)
        self.assertEqual(len(matching.data), 5)

    def test_csr_array(self):
        G, *_ = random_bipartite(5, 5, 0.5, 0)
        matching = maximum_bipartite_matching(G.tocsr())
        self.assertEqual(len(matching.data), 5)

    def test_not_bipartite(self):
        # Complete graph not bipartite, can recognise this by
        # fractional result
        data = [1, 1, 1]
        row = [0, 0, 1]
        col = [1, 2, 2]
        G = sp.coo_array((data, (row, col)), shape=(3, 3))
        with self.assertRaises(ValueError):
            maximum_bipartite_matching(G)


class TestWeightedMatching(unittest.TestCase):
    def test_coo_array(self):
        G, *_ = random_bipartite(5, 5, 20, 0)
        matching = maximum_weighted_matching(G)
        self.assertEqual(len(matching.data), 5)

    def test_csr_array(self):
        G, *_ = random_bipartite(5, 5, 20, 0)
        matching = maximum_weighted_matching(G.tocsr())
        self.assertEqual(len(matching.data), 5)

    def test_not_bipartite(self):
        # Complete graph not bipartite, can recognise this by
        # fractional result
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
