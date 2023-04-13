# Unit tests for your mod. This is separate from the example tests, it should
# test smaller pieces of functionality, and can also test multiple cases.
# This template should be copied to tests/test_<mod-name>.py

import unittest
import numpy as np

# from gurobi_optimods.datasets import load_mod_example_data
from gurobi_optimods.qubo import solve_qubo


class TestQubo(unittest.TestCase):
    # def test_datasets(self):
    #     # If you added something to optimods.datasets, test it here
    #     data = load_mod_example_data()
    #     self.assertEqual(set(data.keys()), {"a", "b", "c"})

    def test_none(self):
        Q = None
        solution = solve_qubo(coeff_matrix=Q)
        self.assertTrue(solution is None)

    def test_positive(self):
        Q = np.array([[2, 5], [3, 6]])
        solution = solve_qubo(coeff_matrix=Q)
        self.assertEqual(solution, np.array([0, 0]))

    def test_negative(self):
        Q = np.array([[-2, -5], [-3, -6]])
        solution = solve_qubo(coeff_matrix=Q)
        self.assertEqual(solution, np.array([1, 1]))
