# Unit tests for your mod. This is separate from the example tests, it should
# test smaller pieces of functionality, and can also test multiple cases.
# This template should be copied to tests/test_<mod-name>.py

import unittest

from gurobi_optimods.datasets import load_mod_example_data
from gurobi_optimods.networkdesign import solve_mod


class TestMod(unittest.TestCase):
    def test_datasets(self):
        # If you added something to optimods.datasets, test it here
        data = load_mod_example_data()
        self.assertEqual(set(data.keys()), {"a", "b", "c"})

    def test_simple(self):
        data = None
        solution = solve_mod(data)
        self.assertTrue(solution is None)
