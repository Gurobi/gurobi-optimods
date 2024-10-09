# Tests that solve models using the public API

import collections
import random
import unittest

import gurobipy as gp

from gurobi_optimods.datasets import load_opf_example
from gurobi_optimods.opf import solve_opf


def size_limited_license():
    with gp.Env(params={"OutputFlag": 0}) as env, gp.Model(env=env) as model:
        model.addVars(2001)
        try:
            model.optimize()
            return False
        except gp.GurobiError:
            return True


class TestInvalidData(unittest.TestCase):
    def setUp(self):
        self.case = load_opf_example("pglib_opf_case5_pjm")

    def test_costtype_1(self):
        # Error out, we only handle costtype 2
        self.case["gencost"][1]["costtype"] = 1
        with self.assertRaisesRegex(
            ValueError, "only generator costtype=2 is supported"
        ):
            solve_opf(self.case, opftype="AC")

    def test_gencost_costvector_length(self):
        # Error out if gencost.n is wrong
        self.case["gencost"][2]["n"] = 4
        with self.assertRaisesRegex(
            ValueError, "mismatch between gencost.n and costvector length"
        ):
            solve_opf(self.case, opftype="AC")

    def test_nonquadratic_cost(self):
        # Error out on any cubic terms or higher
        self.case["gencost"][0]["n"] = 4
        self.case["gencost"][0]["costvector"] = [1, 2, 3, 4]
        with self.assertRaisesRegex(
            ValueError, "only quadratic and linear cost functions"
        ):
            solve_opf(self.case, opftype="AC")

    def test_gencost_matches_gen(self):
        # Error out if not 1:1 between gen and gencost (we don't handle the
        # reactive power case)
        self.case["gen"].append(dict(self.case["gen"][-1]))
        with self.assertRaisesRegex(
            ValueError, "mismatch between gen and gencost records"
        ):
            solve_opf(self.case, opftype="AC")

    def test_bad_branch_fbus(self):
        # All branches must point to valid bus ids
        self.case["branch"][3]["fbus"] = 22
        with self.assertRaisesRegex(
            ValueError, "Unknown bus ID referenced in branch fbus"
        ):
            solve_opf(self.case, opftype="AC")

    def test_bad_branch_tbus(self):
        # All branches must point to valid bus ids
        self.case["branch"][5]["tbus"] = 15
        with self.assertRaisesRegex(
            ValueError, "Unknown bus ID referenced in branch tbus"
        ):
            solve_opf(self.case, opftype="AC")

    def test_bad_gen_bus(self):
        # All generators must point to valid bus ids
        self.case["gen"][1]["bus"] = 48
        with self.assertRaisesRegex(
            ValueError, "Unknown bus ID referenced in generator bus"
        ):
            solve_opf(self.case, opftype="AC")


class TestAPICase5_PJM(unittest.TestCase):
    # Test configurations of case5_pjm with known optimal values. We should
    # get reasonable reproducibility on this small case.

    def setUp(self):
        self.case = load_opf_example("pglib_opf_case5_pjm")

    def assert_approx_equal(self, value, expected, tol):
        self.assertLess(abs(value - expected), tol)

    def assert_solution_valid(self, solution):
        # The solution carries over the same structure as the case input data,
        # and some values (e.g. bus refs) must match exactly in the ordering.

        self.assertEqual(
            set(solution.keys()), set(self.case.keys()) | {"et", "success", "f"}
        )
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

    def test_dc(self):
        solution = solve_opf(self.case, opftype="DC")
        self.assertEqual(solution["success"], 1)
        self.assert_solution_valid(solution)

        self.assert_approx_equal(solution["f"], 17479.89, tol=1e-1)

    def test_ac(self):
        solution = solve_opf(self.case, opftype="AC")
        self.assertEqual(solution["success"], 1)
        self.assert_solution_valid(solution)

        self.assert_approx_equal(solution["f"], 17551.89, tol=1e-1)

    def test_ac_branchswitching(self):
        solution = solve_opf(
            self.case,
            opftype="AC",
            branch_switching=True,
            use_mip_start=False,
            solver_params={"MIPGap": 1e-4},
        )
        self.assertEqual(solution["success"], 1)
        self.assert_solution_valid(solution)

        self.assertIsNotNone(solution["f"])
        self.assert_approx_equal(solution["f"], 15174.06, tol=1e-1)

    def test_ac_relax(self):
        solution = solve_opf(self.case, opftype="ACRelax")
        self.assertEqual(solution["success"], 1)
        self.assert_solution_valid(solution)

        self.assert_approx_equal(solution["f"], 14999.71, tol=1e-1)

    def test_ac_branchswitching_infeasible(self):
        # Modify the case to make it infeasible
        self.case["bus"][1]["Vmax"] = 0.8

        # Solve model, expect failure
        with self.assertRaisesRegex(ValueError, "Infeasible model"):
            solve_opf(self.case, opftype="AC", branch_switching=True)


class TestComputeVoltageAnglesBug(unittest.TestCase):
    # Reordering the buses on input (which does not change the network
    # structure) causes compute_voltage_angles to fail.

    def setUp(self):
        self.case = load_opf_example("pglib_opf_case5_pjm")
        reorder = [2, 1, 3, 4, 0]
        self.case["bus"] = [self.case["bus"][ind] for ind in reorder]

    def assert_approx_equal(self, value, expected, tol):
        self.assertLess(abs(value - expected), tol)

    def test_bug(self):
        # Should run without errors
        solution = solve_opf(self.case, opftype="ac", verbose=False)
        # Only one bus should have zero voltage angle (the reference bus)
        num_zeros = sum(bus["Va"] == 0 for bus in solution["bus"])
        self.assertEqual(num_zeros, 1)


class TestAPICase5_PJMReordered(unittest.TestCase):
    # It should be possible to re-order the input data (i.e. the lists of buses,
    # branches, generators) and get the same result

    def setUp(self):
        self.case = load_opf_example("pglib_opf_case5_pjm")

        self.bus_reorder = [i for i, _ in enumerate(self.case["bus"])]
        self.branch_reorder = [i for i, _ in enumerate(self.case["branch"])]
        self.gen_reorder = [i for i, _ in enumerate(self.case["gen"])]

        random.shuffle(self.bus_reorder)
        random.shuffle(self.branch_reorder)
        random.shuffle(self.gen_reorder)

        self.case_reordered = {
            "baseMVA": self.case["baseMVA"],
            "casename": self.case["casename"],
            "bus": [dict(self.case["bus"][ind]) for ind in self.bus_reorder],
            "branch": [dict(self.case["branch"][ind]) for ind in self.branch_reorder],
            "gen": [dict(self.case["gen"][ind]) for ind in self.gen_reorder],
            "gencost": [dict(self.case["gencost"][ind]) for ind in self.gen_reorder],
        }

    def assert_approx_equal(self, value, expected, tol):
        self.assertLess(abs(value - expected), tol)

    def test_dc(self):
        kwargs = dict(
            opftype="dc",
            solver_params={"OptimalityTol": 1e-6},
        )
        solution_original = solve_opf(self.case, **kwargs)
        solution_reordered = solve_opf(self.case_reordered, **kwargs)

        self.assert_approx_equal(
            solution_original["f"], solution_reordered["f"], tol=1e-1
        )

        for ind, orig_ind in enumerate(self.bus_reorder):
            bus_original = solution_original["bus"][orig_ind]
            bus_reordered = solution_reordered["bus"][ind]
            for key, value in bus_original.items():
                self.assert_approx_equal(bus_reordered[key], value, tol=1e-3)

        for ind, orig_ind in enumerate(self.branch_reorder):
            branch_original = solution_original["branch"][orig_ind]
            branch_reordered = solution_reordered["branch"][ind]
            for key, value in branch_original.items():
                self.assert_approx_equal(branch_reordered[key], value, tol=1e-3)

        for ind, orig_ind in enumerate(self.gen_reorder):
            gen_original = solution_original["gen"][orig_ind]
            gen_reordered = solution_reordered["gen"][ind]
            for key, value in gen_original.items():
                # Some Gen entries are not initialized, e.g., Qc2min, check for nan
                if value == value:
                    self.assert_approx_equal(gen_reordered[key], value, tol=1e-3)

    def test_ac(self):
        # Results are too numerically unstable, just check the call goes through
        # without errors and objective function value is in the ballpark.
        kwargs = dict(
            opftype="ac",
            solver_params={"MIPGap": 1e-4},
        )
        solution_original = solve_opf(self.case, **kwargs)
        solution_reordered = solve_opf(self.case_reordered, **kwargs)

        self.assert_approx_equal(
            solution_original["f"], solution_reordered["f"], tol=1e1
        )

    def test_ac_relax(self):
        # Results are too numerically unstable, just check the call goes through
        # without errors and objective function value is in the ballpark.
        kwargs = dict(
            opftype="acrelax",
            solver_params={"OptimalityTol": 1e-6},
        )
        solution_original = solve_opf(self.case, **kwargs)
        solution_reordered = solve_opf(self.case_reordered, **kwargs)

        self.assert_approx_equal(
            solution_original["f"], solution_reordered["f"], tol=1e1
        )
