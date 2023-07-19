# Tests that solve models using the public API

import collections
import unittest

import gurobipy as gp

from gurobi_optimods.opf import solve_opf_model
from tests.opf import load_case9branchswitching, read_case


def size_limited_license():
    with gp.Env(params={"OutputFlag": 0}) as env, gp.Model(env=env) as model:
        model.addVars(2001)
        try:
            model.optimize()
            return False
        except gp.GurobiError:
            return True


class TestAPICase9(unittest.TestCase):
    # Test configurations of case9 with known optimal values. We should
    # get reasonable reproducibility on this small case.

    def setUp(self):
        self.case = read_case("9")

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
        solution = solve_opf_model(self.case, opftype="DC")
        self.assertEqual(solution["success"], 1)
        self.assert_solution_valid(solution)

        self.assert_approx_equal(solution["f"], 5216.026607, tol=1e-1)
        self.assert_approx_equal(solution["gen"][1]["Pg"], 134.377585, tol=1e-1)
        self.assert_approx_equal(solution["branch"][2]["Pt"], -56.2622, tol=1e-1)

    def test_ac(self):
        solution = solve_opf_model(self.case, opftype="AC")
        self.assertEqual(solution["success"], 1)
        self.assert_solution_valid(solution)

        self.assert_approx_equal(solution["f"], 5296.686204, tol=1e-1)
        self.assert_approx_equal(solution["bus"][2]["Vm"], 1.08662, tol=1e-1)
        self.assert_approx_equal(solution["gen"][1]["Qg"], 0.031844, tol=1e-1)
        self.assert_approx_equal(solution["branch"][0]["Qf"], 12.9656, tol=1e-1)

    def test_ac_branchswitching(self):
        solution = solve_opf_model(
            self.case, opftype="AC", branchswitching=True, usemipstart=False
        )
        self.assertEqual(solution["success"], 1)
        self.assert_solution_valid(solution)

        self.assertIsNotNone(solution["f"])
        self.assert_approx_equal(solution["f"], 5296.6862, tol=1e-1)
        self.assert_approx_equal(solution["bus"][0]["Va"], 0.0, tol=1e-1)
        self.assert_approx_equal(solution["gen"][1]["Qg"], 0.0318441, tol=1e-1)
        self.assert_approx_equal(
            solution["branch"][2]["Pt"], 55.96906046691643, tol=1e-1
        )

    def test_ac_relax(self):
        solution = solve_opf_model(self.case, opftype="ACRelax")
        self.assertEqual(solution["success"], 1)
        self.assert_solution_valid(solution)

        self.assert_approx_equal(solution["f"], 5296.66532, tol=1e-1)
        self.assert_approx_equal(solution["gen"][0]["Pg"], 89.803524, tol=1e-1)
        self.assert_approx_equal(solution["branch"][1]["Pt"], -34.1774, tol=1e-1)

    def test_ac_branchswitching_infeasible(self):
        # Modify the case to make it infeasible
        self.case["bus"][1]["Vmax"] = 0.8

        # Solve model, expect failure
        solution = solve_opf_model(self.case, opftype="AC", branchswitching=True)
        self.assertIsNotNone(solution)
        self.assertEqual(solution["success"], 0)


@unittest.skipIf(size_limited_license(), "size-limited-license")
class TestAPILargeModels(unittest.TestCase):
    # Tests with larger models requiring a full license

    def setUp(self):
        self.cases = ["9", "14", "57", "118", "300"]
        self.casedata = [read_case(case) for case in self.cases]
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
                solution = solve_opf_model(case, opftype="DC")
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
                solution = solve_opf_model(case, opftype="AC")
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
                solution = solve_opf_model(case, opftype="ACRelax")
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


class TestAPIBranchSwitching(unittest.TestCase):
    def setUp(self):
        # Modification of case9 where branch switching has an effect
        self.case = load_case9branchswitching()

    def test_ac_defaults(self):
        # no branch switching, all branches should be on
        for usemipstart in [True, False]:
            with self.subTest(usemipstart=usemipstart):
                solution = solve_opf_model(
                    self.case,
                    opftype="AC",
                    branchswitching=False,
                    usemipstart=usemipstart,
                )
                counts = collections.Counter(
                    branch["switching"] for branch in solution["branch"].values()
                )
                self.assertEqual(counts, {1: 12})

    def test_ac_branchswitching(self):
        # branch switching with no minimum, some branches should be off
        for usemipstart in [True, False]:
            with self.subTest(usemipstart=usemipstart):
                solution = solve_opf_model(
                    self.case,
                    opftype="AC",
                    branchswitching=True,
                    minactivebranches=0.0,
                    usemipstart=usemipstart,
                    solver_params={
                        "MIPGap": 1e-5
                    },  # force to find a solution with branches turned off
                )
                counts = collections.Counter(
                    branch["switching"] for branch in solution["branch"].values()
                )
                self.assertLess(counts[1], 12)
                self.assertGreater(counts[0], 0)

    def test_ac_minactivebranches(self):
        # branch switching with 100% minimum, all branches should be on
        for usemipstart in [True, False]:
            with self.subTest(usemipstart=usemipstart):
                solution = solve_opf_model(
                    self.case,
                    opftype="AC",
                    branchswitching=True,
                    minactivebranches=1.0,
                    usemipstart=usemipstart,
                )
                counts = collections.Counter(
                    branch["switching"] for branch in solution["branch"].values()
                )
                self.assertEqual(counts, {1: 12})

    def test_acrelax_defaults(self):
        # no branch switching, all branches should be on
        for usemipstart in [True, False]:
            with self.subTest(usemipstart=usemipstart):
                solution = solve_opf_model(
                    self.case,
                    opftype="ACRelax",
                    branchswitching=False,
                    usemipstart=usemipstart,
                )
                counts = collections.Counter(
                    branch["switching"] for branch in solution["branch"].values()
                )
                self.assertEqual(counts, {1: 12})

    def test_acrelax_branchswitching(self):
        # branch switching with no minimum, some branches should be off
        for usemipstart in [True, False]:
            with self.subTest(usemipstart=usemipstart):
                solution = solve_opf_model(
                    self.case,
                    opftype="ACRelax",
                    branchswitching=True,
                    minactivebranches=0.0,
                    usemipstart=usemipstart,
                    solver_params={
                        "MIPGap": 0.0
                    },  # need to make sure that we always find those solutions
                )
                counts = collections.Counter(
                    branch["switching"] for branch in solution["branch"].values()
                )
                self.assertLess(counts[1], 12)
                self.assertGreater(counts[0], 0)

    def test_acrelax_minactivebranches(self):
        # branch switching with 100% minimum, all branches should be on
        for usemipstart in [True, False]:
            with self.subTest(usemipstart=usemipstart):
                solution = solve_opf_model(
                    self.case,
                    opftype="ACRelax",
                    branchswitching=True,
                    minactivebranches=1.0,
                    usemipstart=usemipstart,
                )
                counts = collections.Counter(
                    branch["switching"] for branch in solution["branch"].values()
                )
                self.assertEqual(counts, {1: 12})

    def test_dc_defaults(self):
        # no branch switching, all branches should be on
        for usemipstart in [True, False]:
            with self.subTest(usemipstart=usemipstart):
                solution = solve_opf_model(
                    self.case,
                    opftype="DC",
                    branchswitching=False,
                    usemipstart=usemipstart,
                )
                counts = collections.Counter(
                    branch["switching"] for branch in solution["branch"].values()
                )
                self.assertEqual(counts, {1: 12})

    @unittest.expectedFailure
    def test_dc_branchswitching(self):
        # branch switching with no minimum, some branches should be off
        for usemipstart in [True, False]:
            with self.subTest(usemipstart=usemipstart):
                solution = solve_opf_model(
                    self.case,
                    opftype="DC",
                    branchswitching=True,
                    minactivebranches=0.0,
                    usemipstart=usemipstart,
                )
                counts = collections.Counter(
                    branch["switching"] for branch in solution["branch"].values()
                )
                self.assertLess(counts[1], 12)
                self.assertGreater(counts[0], 0)

    def test_dc_minactivebranches(self):
        # branch switching with 100% minimum, all branches should be on
        for usemipstart in [True, False]:
            with self.subTest(usemipstart=usemipstart):
                solution = solve_opf_model(
                    self.case,
                    opftype="DC",
                    branchswitching=True,
                    minactivebranches=0.99,
                    usemipstart=usemipstart,
                )
                counts = collections.Counter(
                    branch["switching"] for branch in solution["branch"].values()
                )
                self.assertEqual(counts, {1: 12})


@unittest.skipIf(size_limited_license(), "size-limited-license")
class TestAPINewYork(unittest.TestCase):
    # test a real data set for New York

    def setUp(self):
        self.case = read_case("NY")

        # New York test values
        self.Va_NY = [
            -4.317424,
            -6.18956,
            -5.959205,
            -6.085699,
            -6.0442,
            -5.796713,
            -5.670924,
            -5.770995,
            -3.930278,
            -5.679304,
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
        self.Pf_NY = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 82.369589]

    def test_dc(self):
        # solve opf model and return a solution
        solution = solve_opf_model(self.case, opftype="DC")
        self.assertIsNotNone(solution)
        self.assertEqual(solution["success"], 1)
        self.assertIsNotNone(solution["f"])
        # differences can be quite big because we solve only to 0.1% optimality
        self.assertLess(abs(solution["f"] - 681590.8014), 1e2)
        for i in range(0, 10):
            self.assertLess(abs(solution["bus"][i + 1]["Va"] - self.Va_NY[i]), 1e1)
            self.assertLess(abs(solution["gen"][i + 1]["Pg"] - self.Pg_NY[i]), 1e1)
            self.assertLess(abs(solution["branch"][i + 1]["Pf"] - self.Pf_NY[i]), 1e1)
