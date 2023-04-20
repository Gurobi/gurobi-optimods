import unittest
import numpy as np
from numpy.testing import assert_array_equal
import scipy.sparse as sp

# from gurobi_optimods.datasets import load_mod_example_data
from gurobi_optimods.qubo import solve_qubo


class TestQubo(unittest.TestCase):
    # def test_datasets(self):
    #     # If you added something to optimods.datasets, test it here
    #     data = load_mod_example_data()
    #     self.assertEqual(set(data.keys()), {"a", "b", "c"})

    def test_none(self):
        Q = None
        val, solution = solve_qubo(coeff_matrix=Q)
        self.assertTrue(val is None)
        self.assertTrue(solution is None)

    def test_nonquadratic(self):
        Q = np.ones((2, 3))
        with self.assertRaises(ValueError):
            solve_qubo(coeff_matrix=Q)

    def test_non2dimensional(self):
        Q = np.ones((1, 3, 4))
        with self.assertRaises(ValueError):
            solve_qubo(coeff_matrix=Q)

    def test_positive(self):
        Q = np.array([[2, 5], [3, 6]])
        val, solution = solve_qubo(coeff_matrix=Q)
        self.assertEqual(val, 0)
        assert_array_equal(solution, np.array([0, 0]))

    def test_negative(self):
        Q = np.array([[-2, -5], [-3, -6]])
        val, solution = solve_qubo(coeff_matrix=Q)
        self.assertEqual(val, -16)
        assert_array_equal(solution, np.array([1, 1]))

    def test_Q1(self):
        Q = np.array([[2, -5], [-3, 4]])
        val, solution = solve_qubo(coeff_matrix=Q)
        self.assertEqual(val, -2)
        assert_array_equal(solution, np.array([1, 1]))

    def test_Q2(self):
        Q = np.array([[1, -2], [2, -3]])
        val, solution = solve_qubo(coeff_matrix=Q)
        self.assertEqual(val, -3)
        assert_array_equal(solution, np.array([0, 1]))
