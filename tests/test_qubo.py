import unittest
import numpy as np
from numpy.testing import assert_array_equal
import scipy.sparse as sp

from gurobi_optimods.qubo import solve_qubo


class TestQubo(unittest.TestCase):
    def test_none(self):
        Q = None
        with self.assertRaises(AttributeError):
            solve_qubo(Q)

    def test_nonquadratic(self):
        Q = np.ones((2, 3))
        with self.assertRaises(ValueError):
            solve_qubo(Q)

    def test_non2dimensional(self):
        Q = np.ones((1, 3, 4))
        with self.assertRaises(ValueError):
            solve_qubo(Q)

    def test_positive(self):
        Q = np.array([[2, 5], [3, 6]])
        result = solve_qubo(Q)
        self.assertEqual(result.objective_value, 0)
        assert_array_equal(result.solution, np.array([0, 0]))

    def test_negative(self):
        Q = np.array([[-2, -5], [-3, -6]])
        result = solve_qubo(Q)
        self.assertEqual(result.objective_value, -16)
        assert_array_equal(result.solution, np.array([1, 1]))

    def test_Q1(self):
        Q = np.array([[2, -5], [-3, 4]])
        result = solve_qubo(Q)
        self.assertEqual(result.objective_value, -2)
        assert_array_equal(result.solution, np.array([1, 1]))

    def test_Q2(self):
        Q = np.array([[1, -2], [2, -3]])
        result = solve_qubo(Q)
        self.assertEqual(result.objective_value, -3)
        assert_array_equal(result.solution, np.array([0, 1]))

    def test_sp(self):
        data = [-1, -2, 3]
        row = [0, 0, 1]
        col = [1, 2, 2]
        Q = sp.coo_matrix((data, (row, col)), shape=(3, 3))
        result = solve_qubo(Q)
        self.assertEqual(result.objective_value, -2)
        assert_array_equal(result.solution, np.array([1, 0, 1]))
