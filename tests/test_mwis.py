import unittest
from itertools import combinations

import numpy as np
import scipy.sparse as sp
from numpy.testing import assert_array_equal

from gurobi_optimods.mwis import (
    maximum_weighted_clique,
    maximum_weighted_independent_set,
)


def get_adjacency_matrix(num_vertices, density, seed, scipy_func=sp.csr_matrix):
    """Create the upper triangular adjacency matrix for a random graph
        with the given number of vertices and density.

    Args:
        num_vertices (int): Number of vertices.
        density (float): Graph density between 0 and 1.
        seed (int): Random seed value.

    Returns:
        matrix: The scipy adjacency matrix/array.
    """
    rng = np.random.default_rng(seed=seed)
    num_edges = int(0.5 * density * num_vertices * (num_vertices - 1))
    if num_edges == 0:
        return scipy_func((num_vertices, num_vertices))
    edges = rng.choice(
        list(combinations(range(num_vertices), 2)), num_edges, replace=False
    )
    rows, cols = zip(*edges)
    data = np.ones(num_edges)
    return scipy_func((data, (rows, cols)), shape=(num_vertices, num_vertices))


class TestInput(unittest.TestCase):
    def test_input_adjacency_matrix_type(self):
        for scipy_func in [
            sp.csr_matrix,
            sp.csc_matrix,
            sp.coo_matrix,
            sp.csr_array,
            sp.csc_array,
            sp.coo_array,
        ]:
            for func in [maximum_weighted_independent_set, maximum_weighted_clique]:
                adjacency_matrix = get_adjacency_matrix(
                    10, 0.75, 0, scipy_func=scipy_func
                )
                func(adjacency_matrix, np.ones(10))

        numpy_matrix = np.triu(np.random.randint(low=0, high=2, size=(10, 10)))
        for func in [maximum_weighted_independent_set, maximum_weighted_clique]:
            with self.assertRaises(ValueError):
                func(numpy_matrix, np.ones(10))

    def test_input_adjacency_upper_triangular(self):
        with_diagonals = get_adjacency_matrix(10, 0.75, 0) + sp.diags(np.ones(10))
        not_triangular = get_adjacency_matrix(10, 0.75, 0) + sp.diags(np.ones(9), -1)
        for adjacency_matrix in [with_diagonals, not_triangular]:
            for func in [maximum_weighted_independent_set, maximum_weighted_clique]:
                with self.assertRaises(ValueError):
                    func(adjacency_matrix, np.ones(10))

    def test_input_weights(self):
        adjacency_matrix = get_adjacency_matrix(10, 0.75, 0)
        weights = [1 for _ in range(10)]
        for func in [maximum_weighted_independent_set, maximum_weighted_clique]:
            with self.assertRaises(ValueError):
                func(adjacency_matrix, weights)


class TestMWIS(unittest.TestCase):
    def test_random_graph(self):
        for density in [np.random.random() for _ in range(5)]:
            with self.subTest(density=density):
                num_vertices, density, seed = 10, density, 0
                adjacency_matrix = get_adjacency_matrix(num_vertices, density, seed)
                weights = np.random.randint(1, 100, size=num_vertices)
                mwis = maximum_weighted_independent_set(adjacency_matrix, weights)
                self.assertLessEqual(len(mwis.x), num_vertices)
                self.assertGreaterEqual(len(mwis.x), 1)
                self.assertLessEqual(mwis.f, weights.sum())

    def test_complete_graph(self):
        num_vertices, density, seed = 10, 1, 0
        adjacency_matrix = get_adjacency_matrix(num_vertices, density, seed)
        weights = np.random.randint(1, 100, size=num_vertices)
        mwis = maximum_weighted_independent_set(adjacency_matrix, weights)
        self.assertEqual(len(mwis.x), 1)
        self.assertEqual(mwis.f, weights.max())

    def test_empty_graph(self):
        num_vertices, density, seed = 10, 0, 0
        adjacency_matrix = get_adjacency_matrix(num_vertices, density, seed)
        weights = np.random.randint(1, 100, size=num_vertices)
        mwis = maximum_weighted_independent_set(adjacency_matrix, weights)
        self.assertEqual(len(mwis.x), num_vertices)
        self.assertEqual(mwis.f, weights.sum())

    def test_known_graph(self):
        rows = [0, 0, 0, 1, 1, 2, 2, 3, 4, 4, 5, 6]
        cols = [1, 3, 4, 3, 5, 3, 6, 7, 5, 7, 6, 7]
        data = [1 for _ in range(12)]
        adjacency_matrix = sp.csr_matrix((data, (rows, cols)), shape=(8, 8))
        weights = np.array([2**i for i in range(8)])
        mwis = maximum_weighted_independent_set(adjacency_matrix, weights)
        assert_array_equal(mwis.x, np.array([0, 2, 5, 7]))
        self.assertEqual(mwis.f, 165)

    def test_solver_params(self):
        num_vertices, density, seed = 10, 0.75, 0
        adjacency_matrix = get_adjacency_matrix(num_vertices, density, seed)
        weights = np.random.randint(1, 100, size=num_vertices)
        mwis = maximum_weighted_independent_set(
            adjacency_matrix, weights, solver_params={"presolve": 2, "method": 1}
        )
        self.assertGreaterEqual(len(mwis.x), 1)
        self.assertLessEqual(mwis.f, weights.sum())


class TestMWC(unittest.TestCase):
    def test_random_graph(self):
        for density in [np.random.random() for _ in range(5)]:
            with self.subTest(density=density):
                num_vertices, density, seed = 10, density, 0
                adjacency_matrix = get_adjacency_matrix(num_vertices, density, seed)
                weights = np.random.randint(1, 100, size=num_vertices)
                mwc = maximum_weighted_clique(adjacency_matrix, weights)
                self.assertLessEqual(len(mwc.x), num_vertices)
                self.assertGreaterEqual(len(mwc.x), 1)
                self.assertLessEqual(mwc.f, weights.sum())

    def test_complete_graph(self):
        num_vertices, density, seed = 10, 1, 0
        adjacency_matrix = get_adjacency_matrix(num_vertices, density, seed)
        weights = np.random.randint(1, 100, size=num_vertices)
        mwc = maximum_weighted_clique(adjacency_matrix, weights)
        self.assertEqual(len(mwc.x), num_vertices)
        self.assertEqual(mwc.f, weights.sum())

    def test_empty_graph(self):
        num_vertices, density, seed = 10, 0, 0
        adjacency_matrix = get_adjacency_matrix(num_vertices, density, seed)
        weights = np.random.randint(1, 100, size=num_vertices)
        mwc = maximum_weighted_clique(adjacency_matrix, weights)
        self.assertEqual(len(mwc.x), 1)
        self.assertEqual(mwc.f, weights.max())

    def test_known_graph(self):
        rows = [0, 0, 0, 1, 1, 2, 2, 3, 4, 4, 5, 6]
        cols = [1, 3, 4, 3, 5, 3, 6, 7, 5, 7, 6, 7]
        data = [1 for _ in range(12)]
        adjacency_matrix = sp.csr_matrix((data, (rows, cols)), shape=(8, 8))
        weights = np.array([2**i for i in range(8)])
        mwc = maximum_weighted_clique(adjacency_matrix, weights)
        assert_array_equal(mwc.x, np.array([6, 7]))
        self.assertEqual(mwc.f, 192)

    def test_solver_params(self):
        num_vertices, density, seed = 10, 0.75, 0
        adjacency_matrix = get_adjacency_matrix(num_vertices, density, seed)
        weights = np.random.randint(1, 100, size=num_vertices)
        mwc = maximum_weighted_independent_set(
            adjacency_matrix, weights, solver_params={"presolve": 2, "method": 1}
        )
        self.assertGreaterEqual(len(mwc.x), 1)
        self.assertLessEqual(mwc.f, weights.sum())


if __name__ == "__main__":
    unittest.main()
