import random
import unittest
from itertools import chain, product

import numpy as np
import scipy.sparse as sp
from numpy.testing import assert_allclose, assert_array_equal

try:
    import networkx as nx
except ImportError:
    nx = None

from gurobi_optimods.matching import (
    maximum_bipartite_matching,
    maximum_weighted_matching,
)


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


class TestBipartiteMatching(unittest.TestCase):
    def assert_is_unweighted_matching(self, matching):
        assert_allclose(matching.data, np.ones(matching.data.shape))
        adj = matching.todense()
        assert_allclose(adj, adj.T)
        self.assertTrue(np.alltrue(adj.sum(axis=0) <= 1))

    def test_empty(self):
        # Matching of an empty graph is empty
        nodes1 = np.arange(5)
        nodes2 = np.arange(5, 13)
        degree = nodes1.shape[0] + nodes2.shape[0]
        adjacency = sp.coo_matrix(([], ([], [])), shape=(degree, degree))

        matching = maximum_bipartite_matching(adjacency, nodes1, nodes2)

        self.assertIsInstance(matching, sp.spmatrix)
        self.assertIsNot(matching, adjacency)
        self.assertEqual(matching.nnz, 0)
        self.assertEqual(matching.shape, (degree, degree))
        self.assert_is_unweighted_matching(matching)

    def test_complete_bipartite(self):
        # Matching of a complete graph has min(n1, n2) edges, and is a matching
        n1, n2 = 6, 8
        nodes1 = np.arange(n1)
        nodes2 = np.arange(n1, n1 + n2)
        degree = nodes1.shape[0] + nodes2.shape[0]
        edges = list(chain(product(nodes1, nodes2), product(nodes2, nodes1)))
        i, j = zip(*edges)
        data = np.ones(len(i))
        adjacency = sp.coo_matrix((data, (i, j)), shape=(degree, degree))
        expect_matching_size = min(n1, n2)

        matching = maximum_bipartite_matching(adjacency, nodes1, nodes2)

        self.assertIsInstance(matching, sp.spmatrix)
        self.assertIsNot(matching, adjacency)
        self.assertEqual(matching.shape, (degree, degree))
        self.assert_is_unweighted_matching(matching)
        self.assertEqual(matching.nnz, expect_matching_size * 2)  # symmetric

    def test_simple(self):
        # Simple case with a known solution
        nodes1 = np.array([1, 3, 4, 6])
        nodes2 = np.array([0, 2, 5])
        degree = nodes1.shape[0] + nodes2.shape[0]
        i = np.array([1, 3, 4, 6])
        j = np.array([2, 5, 5, 2])
        data = np.ones(4)
        triu = sp.coo_matrix((data, (i, j)), shape=(degree, degree))
        adjacency = triu + triu.T
        expect_matching_size = 2

        matching = maximum_bipartite_matching(adjacency, nodes1, nodes2)

        self.assertIsInstance(matching, sp.spmatrix)
        self.assertIsNot(matching, adjacency)
        self.assertEqual(matching.shape, (degree, degree))
        self.assert_is_unweighted_matching(matching)
        self.assertEqual(matching.nnz, expect_matching_size * 2)  # symmetric

    def test_random(self):
        # Property test for matchings on random graphs
        adjacency, nodes1, nodes2 = random_bipartite(n1=2, n2=3, p=0.5, seed=2394)

        matching = maximum_bipartite_matching(adjacency, nodes1, nodes2)

        self.assertIsInstance(matching, sp.spmatrix)
        self.assertIsNot(matching, adjacency)
        self.assertEqual(matching.shape, adjacency.shape)
        self.assert_is_unweighted_matching(matching)

    def test_random_csr_matrix(self):
        # Property test for matchings on random graphs
        adjacency, nodes1, nodes2 = random_bipartite(n1=8, n2=7, p=0.5, seed=98634)

        matching = maximum_bipartite_matching(sp.csr_matrix(adjacency), nodes1, nodes2)

        self.assertIsInstance(matching, sp.spmatrix)
        self.assertIsNot(matching, adjacency)
        self.assertEqual(matching.shape, adjacency.shape)
        self.assert_is_unweighted_matching(matching)

    def test_random_csc_array(self):
        # Property test for matchings on random graphs
        adjacency, nodes1, nodes2 = random_bipartite(n1=5, n2=2, p=0.5, seed=34687)

        matching = maximum_bipartite_matching(sp.csc_array(adjacency), nodes1, nodes2)

        self.assertIsInstance(matching, sp.spmatrix)
        self.assertIsNot(matching, adjacency)
        self.assertEqual(matching.shape, adjacency.shape)
        self.assert_is_unweighted_matching(matching)

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_random_networkx(self):
        graph = nx.bipartite.random_graph(n=5, m=7, p=0.5, seed=293847, directed=True)
        nodes1 = np.array(
            [i for i, node in graph.nodes().items() if node["bipartite"] == 0]
        )
        nodes2 = np.array(
            [i for i, node in graph.nodes().items() if node["bipartite"] == 1]
        )

        matching = maximum_bipartite_matching(graph, nodes1, nodes2)

        self.assertIsInstance(matching, nx.Graph)
        self.assertIsNot(matching, graph)
        self.assertEqual(matching.number_of_nodes(), graph.number_of_nodes())
        self.assert_is_unweighted_matching(
            nx.convert_matrix.to_scipy_sparse_array(matching)
        )


class TestWeightedMatching(unittest.TestCase):
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
