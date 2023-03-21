import unittest

from gurobi_optimods.opf import solve_acopf_model
from gurobi_optimods.datasets import load_opf


class TestOpf(unittest.TestCase):
    def test_simple(self):
        conf, case = load_opf()
        solution = solve_acopf_model(conf, case)
        self.assertTrue(solution is None)
