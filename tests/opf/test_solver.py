# Tests that solve models using the public API

import collections
import math
import os
import random
import sys
import unittest

from gurobipy import GRB

from gurobi_optimods.datasets import load_opf_example
from gurobi_optimods.opf import solve_opf

from ..utils import size_limited_license


class TestInvalidData(unittest.TestCase):
    def setUp(self):
        self.case = load_opf_example("case9")

    def test_costtype_1(self):
        # Error out, we only handle costtype 2
        self.case["gencost"][1]["costtype"] = 1
        with self.assertRaisesRegex(
            ValueError, "only generator costtype=2 is supported"
        ):
            solve_opf(self.case, opftype="ACRGLOBAL")

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


class TestAPICase9(unittest.TestCase):
    # Test configurations of case9 with known optimal values. We should
    # get reasonable reproducibility on this small case.

    def setUp(self):
        self.case = load_opf_example("case9")

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

        self.assert_approx_equal(solution["f"], 5216.026607, tol=1e-1)
        self.assert_approx_equal(solution["gen"][1]["Pg"], 134.377585, tol=1e-1)
        self.assert_approx_equal(solution["branch"][2]["Pt"], -56.2622, tol=1e-1)

    def test_ac(self):
        solution = solve_opf(self.case, opftype="AC")
        self.assertEqual(solution["success"], 1)
        self.assert_solution_valid(solution)

        self.assert_approx_equal(solution["f"], 5296.686204, tol=1e-1)
        self.assert_approx_equal(solution["bus"][2]["Vm"], 1.08662, tol=1e-1)
        self.assert_approx_equal(solution["gen"][1]["Qg"], 0.031844, tol=1e-1)
        self.assert_approx_equal(solution["branch"][0]["Qf"], 12.9656, tol=1e-1)

    @unittest.skip("shaky")
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
        self.assert_approx_equal(solution["f"], 5296.76, tol=1.0)
        self.assert_approx_equal(solution["bus"][0]["Va"], 0.0, tol=1e-1)
        self.assert_approx_equal(solution["gen"][1]["Qg"], 0.0318441, tol=1e-1)
        self.assert_approx_equal(solution["branch"][2]["Pt"], 56.23, tol=1e-1)

    @unittest.skipIf(
        os.environ.get("CI", "false") == "true"
        and sys.platform == "darwin"
        and GRB.VERSION_MAJOR == 12,
        "Numerical trouble on github macos runners",
    )
    def test_ac_relax(self):
        solution = solve_opf(self.case, opftype="ACRelax")
        self.assertEqual(solution["success"], 1)
        self.assert_solution_valid(solution)

        self.assert_approx_equal(solution["f"], 5296.66532, tol=1e-1)
        self.assert_approx_equal(solution["gen"][0]["Pg"], 89.803524, tol=1e-1)
        self.assert_approx_equal(solution["branch"][1]["Pt"], -34.1774, tol=1e-1)

    @unittest.skipIf(GRB.VERSION_MAJOR < 12, "Needs Gurobi 12")
    def test_aclocal(self):
        solution = solve_opf(self.case, opftype="acrlocal")
        self.assertEqual(solution["success"], 1)
        self.assert_solution_valid(solution)

        # A local solution should get somewhere close to optimal
        self.assert_approx_equal(solution["f"], 5296.66532, tol=1e1)

    def test_ac_branchswitching_infeasible(self):
        # Modify the case to make it infeasible
        self.case["bus"][1]["Vmax"] = 0.8

        # Solve model, expect failure
        with self.assertRaisesRegex(ValueError, "Infeasible model"):
            solve_opf(self.case, opftype="ACRGLOBAL", branch_switching=True)


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

    @unittest.skip("shaky")
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
        self.assert_approx_equal(solution["f"], 15174, tol=1)

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
            solve_opf(self.case, opftype="ACRLOCAL", branch_switching=True)

    def test_ac_polar(self):
        kwargs = dict(
            opftype="acpglobal",
            solver_params={"MIPGap": 2e-2, "TimeLimit": 60},
        )
        solution = solve_opf(self.case, **kwargs)
        self.assert_approx_equal(solution["f"], 17551.89, tol=1e-1)

    def test_acp_local(self):
        kwargs = dict(
            opftype="acplocal",
        )
        solution = solve_opf(self.case, **kwargs)
        self.assert_approx_equal(solution["f"], 17551.89, tol=1e-1)

    def test_acr_local(self):
        kwargs = dict(
            opftype="acrlocal",
        )
        solution = solve_opf(self.case, **kwargs)
        self.assert_approx_equal(solution["f"], 17551.89, tol=1e-1)


class TestComputeVoltageAnglesBug(unittest.TestCase):
    # Reordering the buses on input (which does not change the network
    # structure) causes compute_voltage_angles to fail.

    def setUp(self):
        self.case = load_opf_example("case9")
        reorder = [8, 2, 5, 7, 6, 1, 3, 4, 0]
        self.case["bus"] = [self.case["bus"][ind] for ind in reorder]

    def assert_approx_equal(self, value, expected, tol):
        self.assertLess(abs(value - expected), tol)

    def test_bug(self):
        # Should run without errors
        solution = solve_opf(self.case, opftype="acrglobal", verbose=False)
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
                if not math.isnan(value):
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

    @unittest.skipIf(GRB.VERSION_MAJOR < 12, "Needs Gurobi 12")
    def test_aclocal(self):
        # Test AC local option, should arrive at a similar result
        kwargs = dict(opftype="acrlocal")
        solution_original = solve_opf(self.case, **kwargs)
        solution_reordered = solve_opf(self.case_reordered, **kwargs)

        self.assert_approx_equal(
            solution_original["f"], solution_reordered["f"], tol=1e1
        )


@unittest.skipIf(size_limited_license(), "size-limited-license")
class TestAPILargeModels(unittest.TestCase):
    # Tests with larger models requiring a full license

    def setUp(self):
        self.cases = ["case9", "case14", "case57", "case118", "case300"]
        self.casedata = [load_opf_example(case) for case in self.cases]
        # DC test values
        self.objvals_dc = [
            5216.026607,
            7642.591776,
            41006.736942,
            125947.881417,
            706240.290695,
        ]
        self.Va_dc = [6.177764, 6.283185, 6.171413, 5.817455, -5.520424]
        self.Pg_dc = [134.377585, 38.032305, 81.931329, 0, 1.2724979]
        self.Pt_dc = [-56.2622, 69.9608, -18.09715, -102.95381, 24.83]
        # AC test values
        self.objvals_ac = [5296.686204, 8081.187603]
        self.Vm_ac = [1.08662, 1.018801]
        self.Qg_ac = [0.031844, 32.114784]
        self.Qf_ac = [12.9656, -12.67811]
        # AC relaxation test values
        self.objvals_acconv = [
            5296.66532,
            8074.9102,
            41710.3065,
            129338.093,
            718633.78596,
        ]
        self.Pg_acconv = [89.803524, 194.796114, 142.58252, 24.518669, 0.030902]
        self.Pt_acconv = [-34.1774, -71.23414, -29.9637, 23.79936, 56.2152]

    def test_dc(self):
        for i in range(5):
            with self.subTest(case=self.cases[i]):
                case = self.casedata[i]
                solution = solve_opf(case, opftype="DC")
                # Check whether the solution point looks correct. Differences can
                # be quite big because we solve only to 0.1% optimality.
                self.assertIsNotNone(solution)
                self.assertEqual(solution["success"], 1)
                self.assertIsNotNone(solution["f"])
                self.assertLess(abs(solution["f"] - self.objvals_dc[i]), 1e1)
                self.assertLess(abs(solution["bus"][0]["Va"] - self.Va_dc[i]), 1e1)
                self.assertLess(abs(solution["gen"][1]["Pg"] - self.Pg_dc[i]), 1e1)
                self.assertLess(abs(solution["branch"][2]["Pt"] - self.Pt_dc[i]), 1e1)

    def test_ac(self):
        # Exact AC is expensive, so only solve the first two cases.
        for i in range(2):
            with self.subTest(case=self.cases[i]):
                case = self.casedata[i]
                solution = solve_opf(case, opftype="AC")
                # Check whether the solution point looks correct. Differences can
                # be quite big because we solve only to 0.1% optimality.
                self.assertIsNotNone(solution)
                self.assertEqual(solution["success"], 1)
                self.assertIsNotNone(solution["f"])
                self.assertLess(abs(solution["f"] - self.objvals_ac[i]), 1e1)
                self.assertLess(abs(solution["bus"][2]["Vm"] - self.Vm_ac[i]), 1e1)
                self.assertLess(abs(solution["gen"][1]["Qg"] - self.Qg_ac[i]), 1e1)
                self.assertLess(abs(solution["branch"][0]["Qf"] - self.Qf_ac[i]), 1e1)

    def test_ac_relax(self):
        # Case 5 is numerically unstable
        for i in range(4):
            with self.subTest(case=self.cases[i]):
                case = self.casedata[i]
                solution = solve_opf(case, opftype="ACRelax")
                # Check whether the solution point looks correct. Differences can
                # be quite big because bigger models are often numerically unstable.
                self.assertIsNotNone(solution)
                self.assertEqual(solution["success"], 1)
                self.assertIsNotNone(solution["f"])
                self.assertLess(abs(solution["f"] - self.objvals_acconv[i]), 10)
                self.assertLess(abs(solution["gen"][0]["Pg"] - self.Pg_acconv[i]), 10)
                self.assertLess(
                    abs(solution["branch"][1]["Pt"] - self.Pt_acconv[i]), 10
                )


@unittest.skip("unstable")
class TestAPIBranchSwitching(unittest.TestCase):
    def setUp(self):
        # Modification of case9 where branch switching has an effect
        self.case = load_opf_example("case9-switching")

    def test_ac_defaults(self):
        # no branch switching, all branches should be on
        for use_mip_start in [True, False]:
            with self.subTest(use_mip_start=use_mip_start):
                solution = solve_opf(
                    self.case,
                    opftype="AC",
                    branch_switching=False,
                    use_mip_start=use_mip_start,
                )
                counts = collections.Counter(
                    branch["switching"] for branch in solution["branch"]
                )
                self.assertEqual(counts, {1: 12})

    def test_ac_branchswitching(self):
        # branch switching with no minimum, some branches should be off
        for use_mip_start in [True, False]:
            with self.subTest(use_mip_start=use_mip_start):
                solution = solve_opf(
                    self.case,
                    opftype="AC",
                    branch_switching=True,
                    min_active_branches=0.0,
                    use_mip_start=use_mip_start,
                    solver_params={
                        "MIPGap": 1e-5
                    },  # force to find a solution with branches turned off
                )
                counts = collections.Counter(
                    branch["switching"] for branch in solution["branch"]
                )
                self.assertLess(counts[1], 12)
                self.assertGreater(counts[0], 0)

    def test_ac_minactivebranches(self):
        # branch switching with 100% minimum, all branches should be on
        for use_mip_start in [True, False]:
            with self.subTest(use_mip_start=use_mip_start):
                solution = solve_opf(
                    self.case,
                    opftype="AC",
                    branch_switching=True,
                    min_active_branches=1.0,
                    use_mip_start=use_mip_start,
                )
                counts = collections.Counter(
                    branch["switching"] for branch in solution["branch"]
                )
                self.assertEqual(counts, {1: 12})

    def test_acrelax_defaults(self):
        # no branch switching, all branches should be on
        for use_mip_start in [True, False]:
            with self.subTest(use_mip_start=use_mip_start):
                solution = solve_opf(
                    self.case,
                    opftype="ACRelax",
                    branch_switching=False,
                    use_mip_start=use_mip_start,
                )
                counts = collections.Counter(
                    branch["switching"] for branch in solution["branch"]
                )
                self.assertEqual(counts, {1: 12})

    def test_acrelax_branchswitching(self):
        # branch switching with no minimum, some branches should be off
        for use_mip_start in [True, False]:
            with self.subTest(use_mip_start=use_mip_start):
                solution = solve_opf(
                    self.case,
                    opftype="ACRelax",
                    branch_switching=True,
                    min_active_branches=0.0,
                    use_mip_start=use_mip_start,
                    solver_params={
                        "MIPGap": 0.0
                    },  # need to make sure that we always find those solutions
                )
                counts = collections.Counter(
                    branch["switching"] for branch in solution["branch"]
                )
                self.assertLess(counts[1], 12)
                self.assertGreater(counts[0], 0)

    def test_acrelax_minactivebranches(self):
        # branch switching with 100% minimum, all branches should be on
        for use_mip_start in [True, False]:
            with self.subTest(use_mip_start=use_mip_start):
                solution = solve_opf(
                    self.case,
                    opftype="ACRelax",
                    branch_switching=True,
                    min_active_branches=1.0,
                    use_mip_start=use_mip_start,
                )
                counts = collections.Counter(
                    branch["switching"] for branch in solution["branch"]
                )
                self.assertEqual(counts, {1: 12})

    def test_dc_defaults(self):
        # no branch switching, all branches should be on
        for use_mip_start in [True, False]:
            with self.subTest(use_mip_start=use_mip_start):
                solution = solve_opf(
                    self.case,
                    opftype="DC",
                    branch_switching=False,
                    use_mip_start=use_mip_start,
                )
                counts = collections.Counter(
                    branch["switching"] for branch in solution["branch"]
                )
                self.assertEqual(counts, {1: 12})

    @unittest.expectedFailure
    def test_dc_branchswitching(self):
        # branch switching with no minimum, some branches should be off
        for use_mip_start in [True, False]:
            with self.subTest(use_mip_start=use_mip_start):
                solution = solve_opf(
                    self.case,
                    opftype="DC",
                    branch_switching=True,
                    min_active_branches=0.0,
                    use_mip_start=use_mip_start,
                )
                counts = collections.Counter(
                    branch["switching"] for branch in solution["branch"]
                )
                self.assertLess(counts[1], 12)
                self.assertGreater(counts[0], 0)

    def test_dc_minactivebranches(self):
        # branch switching with 100% minimum, all branches should be on
        for use_mip_start in [True, False]:
            with self.subTest(use_mip_start=use_mip_start):
                solution = solve_opf(
                    self.case,
                    opftype="DC",
                    branch_switching=True,
                    min_active_branches=0.99,
                    use_mip_start=use_mip_start,
                )
                counts = collections.Counter(
                    branch["switching"] for branch in solution["branch"]
                )
                self.assertEqual(counts, {1: 12})


@unittest.skipIf(size_limited_license(), "size-limited-license")
class TestAPINewYork(unittest.TestCase):
    # test a real data set for New York

    def setUp(self):
        self.case = load_opf_example("caseNY")

        # New York test values
        self.Va_NY = [
            1.642,
            -0.202,
            0.0122,
            -0.116,
            -0.085,
            0.128,
            0.294,
            0.124,
            1.987,
            0.211,
        ]
        self.Pg_NY = [
            1299.0,
            1012.0,
            45.0,
            45.0,
            1.0,
            882.0,
            655.1,
            1259.3,
            641.8,
            100.0,
        ]
        self.Pf_NY = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 72.0]

    def test_dc(self):
        # solve opf model and return a solution
        solution = solve_opf(
            self.case, opftype="DC", solver_params={"OptimalityTol": 1e-6}
        )
        self.assertIsNotNone(solution)
        self.assertEqual(solution["success"], 1)
        self.assertIsNotNone(solution["f"])
        self.assertLess(abs(solution["f"] - 726388.1), 1e0)
        for i in range(0, 10):
            self.assertLess(abs(solution["bus"][i]["Va"] - self.Va_NY[i]), 1e-1)
            self.assertLess(abs(solution["gen"][i]["Pg"] - self.Pg_NY[i]), 1e-1)
            self.assertLess(abs(solution["branch"][i]["Pf"] - self.Pf_NY[i]), 1e-1)
