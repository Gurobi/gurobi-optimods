import unittest

from itertools import combinations
import scipy.sparse as sp
import numpy as np
from numpy.testing import assert_array_equal

from gurobi_optimods.mwis import maximum_weighted_independent_set


def get_adjacency_matrix(num_vertices, density, seed):
    """Create the upper triangular adjacency matrix for a random graph
        with the given number of vertices and density.

    Args:
        num_vertices (int): Number of vertices.
        density (float): Graph density between 0 and 1.
        seed (int): Random seed value.

    Returns:
        csr_matrix: The adjacency matrix in CSR format.
    """
    rng = np.random.default_rng(seed=seed)
    num_edges = int(0.5 * density * num_vertices * (num_vertices - 1))
    if num_edges == 0:
        return sp.csr_matrix((num_vertices, num_vertices))
    edges = rng.choice(
        list(combinations(range(num_vertices), 2)), num_edges, replace=False
    )
    rows, cols = zip(*edges)
    data = np.ones(num_edges)
    return sp.csr_matrix((data, (rows, cols)), shape=(num_vertices, num_vertices))


class TestMWIS(unittest.TestCase):
    def test_random_graph(self):
        for density in [np.random.random() for _ in range(5)]:
            with self.subTest(density=density):
                num_vertices, density, seed = 10, density, 0
                adjacency_matrix = get_adjacency_matrix(num_vertices, density, seed)
                weights = np.random.randint(1, 100, size=num_vertices)
                mwis = maximum_weighted_independent_set(adjacency_matrix, weights)
                self.assertLessEqual(len(mwis), num_vertices)
                self.assertGreaterEqual(len(mwis), 1)

    def test_complete_graph(self):
        num_vertices, density, seed = 10, 1, 0
        adjacency_matrix = get_adjacency_matrix(num_vertices, density, seed)
        weights = np.random.randint(1, 100, size=num_vertices)
        mwis = maximum_weighted_independent_set(adjacency_matrix, weights)
        self.assertEqual(len(mwis), 1)
        self.assertEqual(sum(weights[mwis]), max(weights))

    def test_empty_graph(self):
        num_vertices, density, seed = 10, 0, 0
        adjacency_matrix = get_adjacency_matrix(num_vertices, density, seed)
        weights = np.random.randint(1, 100, size=num_vertices)
        mwis = maximum_weighted_independent_set(adjacency_matrix, weights)
        self.assertEqual(len(mwis), num_vertices)

    def test_known_graph(self):
        rows = [0, 0, 0, 1, 1, 2, 2, 3, 4, 4, 5, 6]
        cols = [1, 3, 4, 3, 5, 3, 6, 7, 5, 7, 6, 7]
        data = [1 for _ in range(12)]
        adjacency_matrix = sp.csr_matrix((data, (rows, cols)), shape=(8, 8))
        weights = np.array([2**i for i in range(8)])
        mwis = maximum_weighted_independent_set(adjacency_matrix, weights)
        assert_array_equal(mwis, np.array([0, 2, 5, 7]))
        self.assertEqual(sum(weights[mwis]), 165)
