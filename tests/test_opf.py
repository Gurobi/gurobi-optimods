import collections
import unittest

import gurobipy as gp

from gurobi_optimods.datasets import (
    load_case9branchswitching,
    load_case9solution,
    load_caseNYopf,
    load_caseopfmat,
    load_filepath,
    load_opfdictcase,
)
from gurobi_optimods.opf import (
    compute_violations_from_given_voltages,
    read_case_from_mat_file,
    read_coords_from_csv_file,
    read_voltages_from_csv_file,
    solve_opf_model,
)

# If plotly is not installed, some tests will be skipped
try:
    import plotly
except ImportError:
    plotly = None

# If plotly is installed, the opfgraphics module should import ok
if plotly:
    from gurobi_optimods.opf.graphics import (
        generate_opf_solution_figure,
        generate_opf_violations_figure,
    )


def size_limited_license():
    with gp.Env(params={"OutputFlag": 0}) as env, gp.Model(env=env) as model:
        model.addVars(2001)
        try:
            model.optimize()
            return False
        except gp.GurobiError:
            return True


class TestAPIVarious(unittest.TestCase):
    # Various simple cases to test

    def test_infeasible(self):
        case = load_opfdictcase()
        # make case infeasible
        case["bus"][1]["Vmax"] = 0.8
        # solve opf model and return a solution
        solution = solve_opf_model(case, logfile="", opftype="AC", branchswitching=1)
        self.assertIsNotNone(solution)
        self.assertEqual(solution["success"], 0)

    # test reading settings and case file from dicts
    def test_opfdicts(self):
        case = load_opfdictcase()
        solution = solve_opf_model(
            case, opftype="AC", branchswitching=1, usemipstart=False
        )
        # check whether the solution point looks correct
        # differences can be quite big because we solve only to 0.1% optimality
        self.assertIsNotNone(solution)
        self.assertEqual(solution["success"], 1)
        self.assertIsNotNone(solution["f"])
        self.assertLess(abs(solution["f"] - 5296.6862), 1e1)
        self.assertLess(abs(solution["bus"][1]["Va"]), 1)
        self.assertLess(abs(solution["gen"][2]["Qg"] - 0.0318441), 2e1)
        self.assertLess(abs(solution["branch"][3]["Pt"] - 55.96906046691643), 1e1)


class TestAPICase9(unittest.TestCase):
    # Test configurations of case9 with known optimal values

    def setUp(self):
        self.case = read_case_from_mat_file(load_caseopfmat("9"))

    def assert_approx_equal(self, value, expected, tol):
        self.assertLess(abs(value - expected), tol)

    def test_dc(self):
        solution = solve_opf_model(
            self.case,
            opftype="DC",
            branchswitching=0,
            usemipstart=False,
        )
        # Check whether the solution point looks correct. Differences can
        # be quite big because we solve only to 0.1% optimality.
        self.assertIsNotNone(solution)
        self.assertEqual(solution["success"], 1)

        # DC test values
        self.assert_approx_equal(solution["f"], 5216.026607, tol=1e1)
        self.assert_approx_equal(solution["bus"][1]["Va"], 6.177764, tol=1e1)
        self.assert_approx_equal(solution["gen"][2]["Pg"], 134.377585, tol=1e1)
        self.assert_approx_equal(solution["branch"][3]["Pt"], -56.2622, tol=1e1)

    def test_ac(self):
        solution = solve_opf_model(
            self.case,
            opftype="AC",
            polar=False,
            useef=True,
            usejabr=True,
            useactivelossineqs=False,
            branchswitching=0,
            usemipstart=False,
        )
        # Check whether the solution point looks correct. Differences can
        # be quite big because we solve only to 0.1% optimality.
        self.assertIsNotNone(solution)
        self.assertEqual(solution["success"], 1)
        self.assert_approx_equal(solution["f"], 5296.686204, tol=1e1)
        self.assert_approx_equal(solution["bus"][3]["Vm"], 1.08662, tol=1e1)
        self.assert_approx_equal(solution["gen"][2]["Qg"], 0.031844, tol=1e1)
        self.assert_approx_equal(solution["branch"][1]["Qf"], 12.9656, tol=1e1)

    def test_ac_relax(self):
        solution = solve_opf_model(
            self.case,
            opftype="AC",
            polar=False,
            useef=False,
            usejabr=True,
            useactivelossineqs=False,
            branchswitching=0,
            usemipstart=False,
        )
        # Check whether the solution point looks correct. Differences can
        # be quite big because we solve only to 0.1% optimality.
        self.assertIsNotNone(solution)
        self.assertEqual(solution["success"], 1)
        self.assert_approx_equal(solution["f"], 5296.66532, tol=1e1)
        self.assert_approx_equal(solution["gen"][1]["Pg"], 89.803524, tol=1e1)
        self.assert_approx_equal(solution["branch"][2]["Pt"], -34.1774, tol=1e1)


@unittest.skipIf(size_limited_license(), "size-limited-license")
class TestAPILargeModels(unittest.TestCase):
    # Tests with larger models requiring a full license

    def setUp(self):
        self.numcases = 5
        self.cases = ["9", "14", "57", "118", "300"]
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
        for i in range(self.numcases):
            with self.subTest(case=self.cases[i]):
                casefile = load_caseopfmat(self.cases[i])
                case = read_case_from_mat_file(casefile)
                solution = solve_opf_model(
                    case,
                    opftype="DC",
                    branchswitching=0,
                    usemipstart=False,
                )
                # Check whether the solution point looks correct. Differences can
                # be quite big because we solve only to 0.1% optimality.
                self.assertIsNotNone(solution)
                self.assertEqual(solution["success"], 1)
                self.assertIsNotNone(solution["f"])
                self.assertLess(abs(solution["f"] - self.objvals_dc[i]), 1e1)
                self.assertLess(abs(solution["bus"][1]["Va"] - self.Va_dc[i]), 1e1)
                self.assertLess(abs(solution["gen"][2]["Pg"] - self.Pg_dc[i]), 1e1)
                self.assertLess(abs(solution["branch"][3]["Pt"] - self.Pt_dc[i]), 1e1)

    def test_ac(self):
        # Exact AC is expensive, so only solve the first two cases.
        for i in range(2):
            with self.subTest(case=self.cases[i]):
                casefile = load_caseopfmat(self.cases[i])
                case = read_case_from_mat_file(casefile)
                solution = solve_opf_model(
                    case,
                    opftype="AC",
                    polar=False,
                    useef=True,
                    usejabr=True,
                    useactivelossineqs=False,
                    branchswitching=0,
                    usemipstart=False,
                )
                # Check whether the solution point looks correct. Differences can
                # be quite big because we solve only to 0.1% optimality.
                self.assertIsNotNone(solution)
                self.assertEqual(solution["success"], 1)
                self.assertIsNotNone(solution["f"])
                self.assertLess(abs(solution["f"] - self.objvals_ac[i]), 1e1)
                self.assertLess(abs(solution["bus"][3]["Vm"] - self.Vm_ac[i]), 1e1)
                self.assertLess(abs(solution["gen"][2]["Qg"] - self.Qg_ac[i]), 1e1)
                self.assertLess(abs(solution["branch"][1]["Qf"] - self.Qf_ac[i]), 1e1)

    def test_ac_relax(self):
        for i in range(self.numcases):
            with self.subTest(case=self.cases[i]):
                casefile = load_caseopfmat(self.cases[i])
                case = read_case_from_mat_file(casefile)
                solution = solve_opf_model(
                    case,
                    opftype="AC",
                    polar=False,
                    useef=False,
                    usejabr=True,
                    useactivelossineqs=False,
                    branchswitching=0,
                    usemipstart=False,
                )
                # Check whether the solution point looks correct. Differences can
                # be quite big because bigger models are often numerically unstable.
                self.assertIsNotNone(solution)
                self.assertEqual(solution["success"], 1)
                self.assertIsNotNone(solution["f"])
                self.assertLess(abs(solution["f"] - self.objvals_acconv[i]), 10)
                self.assertLess(abs(solution["gen"][1]["Pg"] - self.Pg_acconv[i]), 10)
                self.assertLess(
                    abs(solution["branch"][2]["Pt"] - self.Pt_acconv[i]), 10
                )


class TestAPIBranchSwitching(unittest.TestCase):
    def setUp(self):
        # Modification of case9 where branch switching does something
        self.case = load_case9branchswitching()

    def test_ac_defaults(self):
        # no branch switching, all branches should be on
        for usemipstart in [True, False]:
            with self.subTest(usemipstart=usemipstart):
                solution = solve_opf_model(
                    self.case,
                    opftype="AC",
                    branchswitching=0,
                    polar=False,
                    useef=True,
                    usejabr=True,
                    useactivelossineqs=False,
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
                    branchswitching=1,
                    minactivebranches=0.0,
                    polar=False,
                    useef=True,
                    usejabr=True,
                    useactivelossineqs=False,
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
                    branchswitching=1,
                    minactivebranches=1.0,
                    polar=False,
                    useef=True,
                    usejabr=True,
                    useactivelossineqs=False,
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
                    opftype="AC",
                    branchswitching=0,
                    polar=False,
                    useef=False,
                    usejabr=True,  # JABR is the SOCP relaxation
                    useactivelossineqs=False,
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
                    opftype="AC",
                    branchswitching=1,
                    minactivebranches=0.0,
                    polar=False,
                    useef=False,
                    usejabr=True,
                    useactivelossineqs=False,
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
                    opftype="AC",
                    branchswitching=1,
                    minactivebranches=1.0,
                    polar=False,
                    useef=False,
                    usejabr=False,
                    useactivelossineqs=False,
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
                    branchswitching=0,
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
                    branchswitching=1,
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
                    branchswitching=1,
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
        # load path to case file
        casefile = load_caseNYopf()
        # read case file and return a case dictionary
        case = read_case_from_mat_file(casefile)
        # solve opf model and return a solution
        solution = solve_opf_model(
            case,
            opftype="DC",
        )
        self.assertIsNotNone(solution)
        self.assertEqual(solution["success"], 1)
        self.assertIsNotNone(solution["f"])
        # differences can be quite big because we solve only to 0.1% optimality
        self.assertLess(abs(solution["f"] - 681590.8014), 1e2)
        for i in range(0, 10):
            self.assertLess(abs(solution["bus"][i + 1]["Va"] - self.Va_NY[i]), 1e1)
            self.assertLess(abs(solution["gen"][i + 1]["Pg"] - self.Pg_NY[i]), 1e1)
            self.assertLess(abs(solution["branch"][i + 1]["Pf"] - self.Pf_NY[i]), 1e1)


class TestComputeVoltages(unittest.TestCase):
    # test violation computation out of pre-defined voltage data
    def test_volts(self):
        # get path to csv file holding the input voltages for case 9
        voltsfile = load_filepath("case9volts.csv")
        volts_dict = read_voltages_from_csv_file(voltsfile)
        # load case dictionary
        case = load_opfdictcase()
        # compute violations
        violations = compute_violations_from_given_voltages(case, volts_dict, True)
        self.assertTrue(violations is not None)
        self.assertLess(abs(violations["bus"][2]["Vmviol"]), 1e-9)
        self.assertLess(abs(violations["bus"][4]["Pviol"] + 318.899783), 1e-4)
        self.assertLess(abs(violations["bus"][7]["Qviol"] - 9.8410206), 1e-4)
        self.assertLess(abs(violations["branch"][7]["limitviol"] - 66.33435), 1e-4)


@unittest.skipIf(plotly is None, "plotly is not installed")
class TestGraphicsCase9(unittest.TestCase):
    # Currently, this is just a convenience setting while working on OptiMod
    plot_graphics = False

    def setUp(self):
        # graphics test values
        self.graphics_9_x = [
            1129.2,
            980.2,
            977.6,
            1182.8,
            480.6,
            85.4,
            1079.6,
            528.0,
            0.0,
        ]
        self.graphics_9_y = [
            1066.62,
            132.53,
            220.4,
            0.0,
            777.49,
            569.85,
            1130.71,
            653.08,
            647.86,
        ]

    # test plotting a solution from pre-loaded data
    def test_graphics(self):
        # get path to csv file holding the coordinates for case 9
        coordsfile = load_filepath("case9coords.csv")
        coords_dict = read_coords_from_csv_file(coordsfile)
        # load case dictionary
        case = load_opfdictcase()
        # load a precomputed solution and objective value
        solution = load_case9solution()
        # plot the given solution
        fig = generate_opf_solution_figure(case, coords_dict, solution)
        # check whether figure coordinates and scaled input coordinates are the same
        for i in range(9):
            self.assertLess(abs(fig.data[1].x[i] - self.graphics_9_x[i]), 1e-9)
            self.assertLess(abs(fig.data[1].y[i] - self.graphics_9_y[i]), 1e-9)
        if self.plot_graphics:
            fig.show()

    # test plotting a solution from pre-loaded data
    def test_graphics_volts(self):
        # get path to csv file holding the voltage information for case 9
        voltsfile = load_filepath("case9volts.csv")
        volts_dict = read_voltages_from_csv_file(voltsfile)
        # load case dictionary
        case = load_opfdictcase()
        # compute violations
        violations = compute_violations_from_given_voltages(case, volts_dict, True)
        self.assertTrue(violations is not None)

        coordsfile = load_filepath("case9coords.csv")
        coords_dict = read_coords_from_csv_file(coordsfile)
        fig = generate_opf_violations_figure(case, coords_dict, violations)
        if self.plot_graphics:
            fig.show()

    # test plotting a solution after optimization is performed
    def test_graphics_after_solving(self):
        # load case dictionary
        case = load_opfdictcase()
        # solve opf model and return a solution
        solution = solve_opf_model(case, opftype="AC")
        # plot the computed solution
        coordsfile = load_filepath("case9coords.csv")
        coords_dict = read_coords_from_csv_file(coordsfile)
        fig = generate_opf_solution_figure(case, coords_dict, solution)
        # check whether figure coordinates and scaled input coordinates are the same
        for i in range(9):
            self.assertLess(abs(fig.data[1].x[i] - self.graphics_9_x[i]), 1e-9)
            self.assertLess(abs(fig.data[1].y[i] - self.graphics_9_y[i]), 1e-9)
        if self.plot_graphics:
            fig.show()

    # test plotting a solution from pre-loaded data
    def test_graphics_branchswitching(self):
        # get path to csv file holding the coordinates for case 9
        coordsfile = load_filepath("case9coords.csv")
        coords_dict = read_coords_from_csv_file(coordsfile)
        # load case dictionary
        case = load_case9branchswitching()
        # compute a solution
        solution = solve_opf_model(case, opftype="AC", branchswitching=1)
        counts = collections.Counter(
            branch["switching"] for branch in solution["branch"].values()
        )
        self.assertEqual(counts[1], 10)
        self.assertEqual(counts[0], 2)
        # plot the given solution
        fig = generate_opf_solution_figure(case, coords_dict, solution)
        if self.plot_graphics:
            fig.show()


@unittest.skipIf(plotly is None, "plotly is not installed")
@unittest.skipIf(size_limited_license(), "size-limited-license")
class TestGraphicsNewYork(unittest.TestCase):
    # Currently, this is just a convenience setting while working on OptiMod
    plot_graphics = False

    # test a real data set for New York
    def test_NY_graphics(self):
        # load path to case file
        casefile = load_caseNYopf()
        # read case file and return a case dictionary
        case = read_case_from_mat_file(casefile)
        # solve opf model and return a solution
        solution = solve_opf_model(case, opftype="DC")
        self.assertIsNotNone(solution)
        self.assertEqual(solution["success"], 1)
        self.assertIsNotNone(solution["f"])

        if self.plot_graphics:
            # get path to csv file holding the coordinates for NY
            coordsfile = load_filepath("nybuses.csv")
            coords_dict = read_coords_from_csv_file(coordsfile)
            # plot the given solution
            fig = generate_opf_solution_figure(case, coords_dict, solution)
            # test a few coordinates
            self.assertLess(abs(fig.data[1].x[0] - 1381.2), 1e-9)
            self.assertLess(abs(fig.data[1].y[0] - 1203.5), 1e-9)
            self.assertLess(abs(fig.data[1].x[-1] - 837.2), 1e-9)
            self.assertLess(abs(fig.data[1].y[-1] - 511.85), 1e-9)
            fig.show()

    def test_NY_branchswitching(self):
        # load path to case file
        casefile = load_caseNYopf()
        # read case file and return a case dictionary
        case = read_case_from_mat_file(casefile)

        # solve opf model and return a solution
        solution = solve_opf_model(
            case,
            opftype="DC",
            branchswitching=True,
            minactivebranches=0.999,
            solver_params={"TimeLimit": 1},
        )
        self.assertIsNotNone(solution)
        self.assertEqual(solution["success"], 1)
        self.assertIsNotNone(solution["f"])

        # get path to csv file holding the coordinates for NY
        coordsfile = load_filepath("nybuses.csv")
        coords_dict = read_coords_from_csv_file(coordsfile)
        fig = generate_opf_solution_figure(case, coords_dict, solution)
        if self.plot_graphics:
            fig.show()


@unittest.skip("Tests internal options; not needed in CI")
class TestInternal(unittest.TestCase):
    # Test internal options we haven't exposed on the public API yet

    def test_ac_polar(self):
        # Test polar formulation. It's expensive to test all combinations,
        # so for now just test the defaults.
        casefile = load_caseopfmat("9")
        case = read_case_from_mat_file(casefile)
        solution = solve_opf_model(
            case,
            opftype="AC",
            polar=True,
            solver_params={"SolutionLimit": 1},
        )
        self.assertIsNotNone(solution)
        self.assertEqual(solution["success"], 1)

    def test_ac_settings(self):
        # Test all AC solver options with polar=False.
        settingslist = [
            dict(
                opftype="AC",
                polar=False,
                useef=useef,
                usejabr=usejabr,
                branchswitching=branchswitching,
                usemipstart=usemipstart,
                useactivelossineqs=useactivelossineqs,
                minactivebranches=minactivebranches,
            )
            for useef in [False, True]
            for usejabr in [False, True]
            for branchswitching in [0, 1, 2]
            for usemipstart in [False, True]
            for useactivelossineqs in [False, True]
            for minactivebranches in [0, 0.5, 0.95]
        ]

        # Solve opf model for the same case and return a solution.
        # The model has to be feasible for every setting combination.
        casefile = load_caseopfmat("9")
        case = read_case_from_mat_file(casefile)
        for acopf_settings in settingslist:
            with self.subTest(acopf_settings):
                solution = solve_opf_model(
                    case,
                    solver_params={"SolutionLimit": 1},
                    **acopf_settings,
                )
                self.assertIsNotNone(solution)
                self.assertEqual(solution["success"], 1)

    def test_dc_settings(self):
        # DC has more limited configurations: polar, ef, jabr,
        # branchswitching=2, and ivtype are not available
        settingslist = [
            dict(
                opftype="DC",
                branchswitching=branchswitching,
                usemipstart=usemipstart,
                useactivelossineqs=useactivelossineqs,
                minactivebranches=minactivebranches,
            )
            for branchswitching in [0, 1]
            for usemipstart in [False, True]
            for useactivelossineqs in [False, True]
            for minactivebranches in [0, 0.5, 0.95]
        ]

        # Solve opf model for the same case and return a solution.
        # The model has to be feasible for every setting combination.
        casefile = load_caseopfmat("9")
        case = read_case_from_mat_file(casefile)
        for acopf_settings in settingslist:
            with self.subTest(acopf_settings):
                solution = solve_opf_model(
                    case,
                    solver_params={"SolutionLimit": 1},
                    **acopf_settings,
                )
                self.assertIsNotNone(solution)
                self.assertEqual(solution["success"], 1)

    def test_iv_settings(self):
        # IV has more limited configurations: polar, ef, jabr,
        # branchswitching, and mipstart are not available
        settingslist = [
            dict(
                opftype="IV",
                ivtype=ivtype,
                useactivelossineqs=useactivelossineqs,
                minactivebranches=minactivebranches,
            )
            for ivtype in ["plain", "aggressive"]
            for useactivelossineqs in [False, True]
            for minactivebranches in [0, 0.5, 0.95]
        ]

        # Solve opf model for the same case and return a solution.
        # The model has to be feasible for every setting combination.
        casefile = load_caseopfmat("9")
        case = read_case_from_mat_file(casefile)
        for acopf_settings in settingslist:
            with self.subTest(acopf_settings):
                solution = solve_opf_model(
                    case,
                    solver_params={"SolutionLimit": 1},
                    **acopf_settings,
                )
                self.assertIsNotNone(solution)
                self.assertEqual(solution["success"], 1)

    # test IV formulation
    def test_ivopf(self):
        # currently all other cases take very long in IV formulation
        casefile = load_caseopfmat("9")
        # read case file and return a case dictionary
        case = read_case_from_mat_file(casefile)
        # solve opf model and return a solution
        solution = solve_opf_model(case, opftype="IV", ivtype="aggressive")
        # check whether the solution points looks correct
        self.assertIsNotNone(solution)
        self.assertEqual(solution["success"], 1)
        self.assertIsNotNone(solution["f"])
        # differences can be quite big because we solve only to 0.1% optimality
        self.assertLess(abs(solution["f"] - 5297.0142014), 1e1)
        self.assertLess(abs(solution["bus"][1]["Vm"] - 1.0986917), 1e-1)
        self.assertLess(abs(solution["gen"][2]["Pg"] - 132.87813), 1e1)
        self.assertLess(abs(solution["gen"][3]["Qg"] + 22.802347), 1e1)
        self.assertLess(abs(solution["branch"][4]["Pf"] - 95.113306), 1e1)
        self.assertLess(abs(solution["branch"][5]["Qt"] + 18.431373), 1e1)
