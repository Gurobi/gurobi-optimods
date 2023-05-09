import unittest

import gurobipy as gp

from gurobi_optimods.opf import (
    solve_opf_model,
    read_coords_from_csv_file,
    read_voltages_from_csv_file,
    read_case_from_mat_file,
    compute_violations_from_given_voltages,
)
from gurobi_optimods.datasets import (
    load_caseopfmat,
    load_caseNYopf,
    load_opfdictcase,
    load_filepath,
    load_case9solution,
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


class TestOpf(unittest.TestCase):

    numcases = 5
    cases = ["9", "14", "57", "118", "300"]
    # DC test values
    objvals_dc = [5216.026607, 7642.591776, 41006.736942, 125947.881417, 706240.290695]
    Va_dc = [6.177764, 6.283185, 6.171413, 5.817455, -5.520424]
    Pg_dc = [134.377585, 38.032305, 81.931329, 0, 1.2724979]
    Pt_dc = [-56.2622, 69.9608, -18.09715, -102.95381, 24.83]
    # AC test values
    objvals_ac = [5296.686204, 8081.187603]
    Vm_ac = [1.08662, 1.018801]
    Qg_ac = [0.031844, 32.114784]
    Qf_ac = [12.9656, -12.67811]
    # AC relaxation test values
    objvals_acconv = [5296.66532, 8074.9102, 41710.3065, 129338.093, 718633.78596]
    Pg_acconv = [89.803524, 194.796114, 142.58252, 24.518669, 0.030902]
    Pt_acconv = [-34.1774, -71.23414, -29.9637, 23.79936, 56.2152]
    # New York test values
    Va_NY = [
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
    Pg_NY = [1299.0, 1012.0, 45.0, 45.0, 1.0, 882.0, 655.1, 1259.3, 641.8, 100.0]
    Pf_NY = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 82.369589]

    # test simple is on purpose the same as test_acopf for now
    # will be removed in final version
    def test_simple(self):
        # load path to case file
        casefile = load_caseopfmat("9")
        # read case file and return a case dictionary
        case = read_case_from_mat_file(casefile)
        # solve opf model and return a solution
        solution = solve_opf_model(
            case, logfile="", opftype="AC", useef=True, usejabr=True, branchswitching=1
        )
        # check whether the solution points looks correct
        self.assertTrue(solution is not None)
        self.assertTrue(solution["success"] == 1)
        self.assertTrue(solution["f"] is not None)

    # test all possible combinations of user-relevant settings
    @unittest.skip(
        "Skipping test_settings, because it takes too much time. It should be run manually"
    )
    def test_settings(self):
        # construct all possible combinations of user-relevant settings
        settingslist = [
            (
                type,
                polar,
                ef,
                usejabr,
                ivtype,
                branchswitching,
                usemipstart,
                useactivelossineqs,
                minactivebranches,
            )
            for type in ["AC", "DC", "IV"]
            for polar in [False, True]
            for ef in [False, True]
            for usejabr in [False, True]
            for ivtype in ["plain", "aggressive"]
            for branchswitching in [0, 1, 2]
            for usemipstart in [False, True]
            for useactivelossineqs in [False, True]
            for minactivebranches in [0, 0.5, 0.95]
        ]
        # load path to case file
        casefile = load_caseopfmat("9")
        # read case file and return a case dictionary
        case = read_case_from_mat_file(casefile)
        # solve opf model and return a solution
        for s in settingslist:
            # gurobi param.prm file should read
            # SolutionLimit 1
            # MIPGap 0.01
            solution = solve_opf_model(
                case,
                logfile="",
                opftype=s[0],
                polar=s[1],
                useef=s[2],
                usejabr=s[3],
                ivtype=s[4],
                branchswitching=s[5],
                usemipstart=s[6],
                useactivelossineqs=s[7],
                minactivebranches=s[8],
                additional_settings={"gurobiparamfile": "param.prm"},
            )
            print(
                f"running setting opftype={s[0]}, polar={s[1]}, useef={s[2]}, usejabr={s[3]}, ivtype={s[4]},\n"
                f"branchswitching={s[5]}, usemipstart={s[6]}, useactivelossineqs={s[7]}, minactivebranches={s[8]}"
            )
            print(s)
            self.assertTrue(solution is not None)
            self.assertTrue(
                solution["success"] == 1
            )  # the model has to be feasible for every setting combination

    def test_infeasible(self):
        case = load_opfdictcase()
        # make case infeasible
        case["bus"][1]["Vmax"] = 0.8
        # solve opf model and return a solution
        solution = solve_opf_model(case, logfile="", opftype="AC", branchswitching=1)
        self.assertTrue(solution is not None)
        self.assertTrue(solution["success"] == 0)

    # test a real data set for New York
    @unittest.skipIf(size_limited_license(), "size-limited-license")
    def test_NY(self):
        # load path to case file
        casefile = load_caseNYopf()
        # read case file and return a case dictionary
        case = read_case_from_mat_file(casefile)
        # solve opf model and return a solution
        solution = solve_opf_model(case, opftype="DC")
        self.assertTrue(solution is not None)
        self.assertTrue(solution["success"] == 1)
        self.assertTrue(solution["f"] is not None)
        # differences can be quite big because we solve only to 1% optimality
        self.assertLess(abs(solution["f"] - 681590.8014), 1e2)
        for i in range(0, 10):
            self.assertLess(abs(solution["bus"][i + 1]["Va"] - self.Va_NY[i]), 1e1)
            self.assertLess(abs(solution["gen"][i + 1]["Pg"] - self.Pg_NY[i]), 1e1)
            self.assertLess(abs(solution["branch"][i + 1]["Pf"] - self.Pf_NY[i]), 1e1)

    # test reading settings and case file from dicts
    def test_opfdicts(self):
        case = load_opfdictcase()
        solution = solve_opf_model(
            case, opftype="AC", branchswitching=1, usemipstart=False
        )
        # check whether the solution point looks correct
        # differences can be quite big because we solve only to 1% optimality
        self.assertTrue(solution is not None)
        self.assertTrue(solution["success"] == 1)
        self.assertTrue(solution["f"] is not None)
        self.assertLess(abs(solution["f"] - 5296.6862), 1e1)
        self.assertLess(abs(solution["bus"][1]["Va"]), 1)
        self.assertLess(abs(solution["gen"][2]["Qg"] - 0.0318441), 2e1)
        self.assertLess(abs(solution["branch"][3]["Pt"] - 55.96906046691643), 1e1)

    # test DC formulation
    @unittest.skipIf(size_limited_license(), "size-limited-license")
    def test_dcopf(self):

        for i in range(self.numcases):
            # load path to case file in .m and .mat format
            casefile = load_caseopfmat(self.cases[i])
            # read case file and .mat format and return a case dictionary
            case = read_case_from_mat_file(casefile)
            # solve opf models and return a solution
            solution = solve_opf_model(case, opftype="DC", branchswitching=1)
            # check whether the solution point looks correct
            self.assertTrue(solution is not None)
            self.assertTrue(solution["success"] == 1)
            self.assertTrue(solution["f"] is not None)
            # differences can be quite big because we solve only to 1% optimality
            self.assertLess(abs(solution["f"] - self.objvals_dc[i]), 1e1)
            self.assertLess(abs(solution["bus"][1]["Va"] - self.Va_dc[i]), 1e1)
            self.assertLess(abs(solution["gen"][2]["Pg"] - self.Pg_dc[i]), 1e1)
            self.assertLess(abs(solution["branch"][3]["Pt"] - self.Pt_dc[i]), 1e1)

    # test AC formulation
    @unittest.skipIf(size_limited_license(), "size-limited-license")
    def test_acopf(self):

        for i in range(2):
            # load path to case file in .m and .mat format
            casefile = load_caseopfmat(self.cases[i])
            # read case file and .mat format and return a case dictionary
            case = read_case_from_mat_file(casefile)
            # solve opf models and return a solution
            solution = solve_opf_model(case, opftype="Ac")
            # check whether the solution point looks correct
            self.assertTrue(solution is not None)
            self.assertTrue(solution["success"] == 1)
            self.assertTrue(solution["f"] is not None)
            # differences can be quite big because we solve only to 1% optimality
            self.assertLess(abs(solution["f"] - self.objvals_ac[i]), 1e1)
            self.assertLess(abs(solution["bus"][3]["Vm"] - self.Vm_ac[i]), 1e1)
            self.assertLess(abs(solution["gen"][2]["Qg"] - self.Qg_ac[i]), 1e1)
            self.assertLess(abs(solution["branch"][1]["Qf"] - self.Qf_ac[i]), 1e1)

    # test AC formulation relaxation
    @unittest.skipIf(size_limited_license(), "size-limited-license")
    def test_acopfconvex(self):

        for i in range(self.numcases):
            # load path to case file in .m and .mat format
            casefile = load_caseopfmat(self.cases[i])
            # read case file in .mat format and return a case dictionary
            case = read_case_from_mat_file(casefile)
            # solve opf models and return a solution
            solution = solve_opf_model(case, opftype="AC", useef=False)
            # check whether the solution point looks correct
            # differences can be quite big because bigger models are often numerically unstable
            self.assertTrue(solution is not None)
            self.assertTrue(solution["success"] == 1)
            self.assertTrue(solution["f"] is not None)
            self.assertLess(abs(solution["f"] - self.objvals_acconv[i]), 10)
            self.assertLess(abs(solution["gen"][1]["Pg"] - self.Pg_acconv[i]), 10)
            self.assertLess(abs(solution["branch"][2]["Pt"] - self.Pt_acconv[i]), 10)

    # test IV formulation
    def test_ivopf(self):
        # currently all other cases take very long in IV formulation
        casefile = load_caseopfmat("9")
        # read case file and return a case dictionary
        case = read_case_from_mat_file(casefile)
        # solve opf model and return a solution
        solution = solve_opf_model(case, opftype="IV", ivtype="aggressive")
        # check whether the solution points looks correct
        self.assertTrue(solution is not None)
        self.assertTrue(solution["success"] == 1)
        self.assertTrue(solution["f"] is not None)
        # differences can be quite big because we solve only to 1% optimality
        self.assertLess(abs(solution["f"] - 5297.0142014), 1e1)
        self.assertLess(abs(solution["gen"][2]["Pg"] - 132.87813), 1e1)
        self.assertLess(abs(solution["gen"][3]["Qg"] + 22.802347), 1e1)
        self.assertLess(abs(solution["branch"][4]["Pf"] - 95.113306), 1e1)
        self.assertLess(abs(solution["branch"][5]["Qt"] + 18.431373), 1e1)

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
class TestOpfGraphics(unittest.TestCase):

    # Currently, this is just a convenience setting while working on OptiMod
    plot_graphics = False

    # graphics test values
    graphics_9_x = [1129.2, 980.2, 977.6, 1182.8, 480.6, 85.4, 1079.6, 528.0, 0.0]
    graphics_9_y = [
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
        # get path to csv file holding the coordinates for case 9
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
        solution = solve_opf_model(case, opftype="AC", branchswitching=1)
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

    # test a real data set for New York
    @unittest.skipIf(size_limited_license(), "size-limited-license")
    def test_NY_graphics(self):
        # load path to case file
        casefile = load_caseNYopf()
        # read case file and return a case dictionary
        case = read_case_from_mat_file(casefile)
        # solve opf model and return a solution
        solution = solve_opf_model(case, opftype="DC")
        self.assertTrue(solution is not None)
        self.assertTrue(solution["success"] == 1)
        self.assertTrue(solution["f"] is not None)

        self.plot_graphics = True
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

    @unittest.skip(
        "Skipping test_NY_branchswitching, because it takes too much time. It should be run manually.\
        Or we have to pass gurobi settings with a TimeLimit via a .prm file."
    )
    def test_NY_branchswitching(self):
        # load path to case file
        casefile = load_caseNYopf()
        # read case file and return a case dictionary
        case = read_case_from_mat_file(casefile)

        # solve opf model and return a solution
        solution = solve_opf_model(
            case, opftype="DC", branchswitching=True, minactivebranches=0.999
        )
        self.assertTrue(solution is not None)
        self.assertTrue(solution["success"] == 1)
        self.assertTrue(solution["f"] is not None)

        # get path to csv file holding the coordinates for NY
        coordsfile = load_filepath("nybuses.csv")
        coords_dict = read_coords_from_csv_file(coordsfile)
        fig = generate_opf_solution_figure(case, coords_dict, solution)
        if self.plot_graphics:
            fig.show()
