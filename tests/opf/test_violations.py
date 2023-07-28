# Tests of the voltage/violations public API functions

import unittest

from gurobi_optimods.datasets import load_opf_example, load_opf_extra
from gurobi_optimods.opf import compute_violations


class TestComputeVoltages(unittest.TestCase):
    # test violation computation out of pre-defined voltage data

    def setUp(self):
        self.case = load_opf_example("case9")
        self.volts_data = load_opf_extra("case9-voltages")

    def assert_approx_equal(self, value, expected, tol):
        self.assertLess(abs(value - expected), tol)

    def assert_solution_valid(self, solution):
        # The solution carries over the same structure as the case input data,
        # and some values (e.g. bus refs) must match exactly in the ordering.

        self.assertEqual(set(solution.keys()), set(self.case.keys()))
        self.assertEqual(solution["baseMVA"], self.case["baseMVA"])

        for sol_bus, case_bus in zip(solution["bus"], self.case["bus"]):
            self.assertTrue(set(sol_bus.keys()).issuperset((set(case_bus.keys()))))
            self.assertEqual(sol_bus["bus_i"], case_bus["bus_i"])

        for sol_branch, case_branch in zip(solution["branch"], self.case["branch"]):
            self.assertTrue(
                set(sol_branch.keys()).issuperset((set(case_branch.keys())))
            )
            self.assertEqual(sol_branch["fbus"], case_branch["fbus"])
            self.assertEqual(sol_branch["tbus"], case_branch["tbus"])

        for sol_gen, case_gen in zip(solution["gen"], self.case["gen"]):
            self.assertTrue(set(sol_gen.keys()).issuperset((set(case_gen.keys()))))
            self.assertEqual(sol_gen["bus"], case_gen["bus"])

        for sol_gencost, case_gencost in zip(solution["gencost"], self.case["gencost"]):
            self.assertTrue(
                set(sol_gencost.keys()).issuperset((set(case_gencost.keys())))
            )
            self.assertEqual(sol_gencost["costtype"], case_gencost["costtype"])
            self.assertEqual(sol_gencost["n"], case_gencost["n"])

    def test_volts(self):
        violations = compute_violations(self.case, self.volts_data, polar=True)
        self.assert_solution_valid(violations)

        # Known values, should be stable
        self.assert_approx_equal(violations["bus"][1]["Vmviol"], 0.0, tol=1e-9)
        self.assert_approx_equal(violations["bus"][3]["Pviol"], -318.899783, tol=1e-4)
        self.assert_approx_equal(violations["bus"][6]["Qviol"], 9.8410206, tol=1e-4)
        self.assert_approx_equal(
            violations["branch"][6]["limitviol"], 66.33435, tol=1e-4
        )
