import unittest

from gurobi_optimods.opf import (
    solve_opf_model,
    plot_opf_solution,
    read_settings_from_file,
    read_case_from_file,
    read_case_from_mat_file,
)
from gurobi_optimods.datasets import (
    load_caseopf,
    load_caseopfmat,
    load_caseNYopf,
    load_acopfsettings,
    load_dcopfsettings,
    load_ivopfsettings,
    load_simpleopfsettings,
    load_opfdictcase,
    load_opfdictsettings,
    load_opfdictgraphicssettings,
    load_opfgraphicssettings,
    load_case9solution,
)


class TestOpf(unittest.TestCase):

    numcases = 5
    cases = ["9", "14", "57", "118", "300"]
    numvars_dc = [50, 95, 363, 850, 1904]
    objvals_dc = [5216.026607, 7642.591776, 41006.736942, 125947.881417, 706240.290695]
    lincost_dc = [688.133508, 5179.999999, 29931.879752, 84840, 470847.849469]
    theta_9_dc = [6.083, 6.008502, 6.177336, 6.12929, -5.579928]

    # test simple is on purpose the same as test_acopf for now
    # will be removed in final version
    def test_simple(self):
        # load path to settings file
        settingsfile = load_simpleopfsettings()
        # read settings file and return a settings dictionary
        settings = read_settings_from_file(settingsfile)
        # load path to case file
        casefile = load_caseopf("9")
        # read case file and return a case dictionary
        case = read_case_from_file(casefile)
        # solve opf model and return a solution and the final objective value
        solution, objval = solve_opf_model(settings, case, "OPF.log")
        # check whether the solution points looks correct
        self.assertTrue(solution is not None)
        self.assertTrue(objval is not None)

    # test a real data set for New York
    def test_NY(self):
        settings = {"dodc": True}
        # load path to case file
        casefile, casemat = load_caseNYopf()
        # read case file and return a case dictionary
        case = read_case_from_file(casefile)
        # solve opf model and return a solution and the final objective value
        solution1, objval1 = solve_opf_model(settings, case, "OPF.log")
        self.assertTrue(solution1 is not None)
        self.assertTrue(objval1 is not None)

        # read mat file and return a case dictionary
        case = read_case_from_mat_file(casemat)
        # solve opf model and return a solution and the final objective value
        solution2, objval2 = solve_opf_model(settings, case, "OPF.log")
        self.assertTrue(solution2 is not None)
        self.assertTrue(objval2 is not None)

        # solutions should be the same because it's the same data
        self.assertTrue(objval1 == objval2)
        for s in solution1:
            self.assertTrue(solution1[s] == solution2[s])

    # test reading settings and case file from dicts
    def test_opfdicts(self):
        settings = load_opfdictsettings()
        case = load_opfdictcase()
        solution, objval = solve_opf_model(settings, case)
        # check whether the solution points looks correct
        self.assertTrue(solution is not None)
        self.assertTrue(objval is not None)
        self.assertTrue(len(solution) == 134)
        self.assertLess(abs(objval - 5296.665647261), 1e-4)
        self.assertTrue("lincost" in solution.keys())
        self.assertLess(abs(solution["lincost"] - 704.399018), 1e-4)
        self.assertTrue("twinP_1_1_4" in solution.keys())
        self.assertLess(abs(solution["twinP_1_1_4"]), 1e-4)
        self.assertTrue("Q_9_4_9" in solution.keys())
        self.assertLess(abs(solution["Q_9_4_9"] - 0.127732), 1e-4)

    # test mat file
    def test_matfile(self):
        # load path to settings file
        settingsfile = load_dcopfsettings()
        # read settings file and return a settings dictionary
        settings = read_settings_from_file(settingsfile)
        # load path to case file

        for i in range(self.numcases):
            # load path to case file
            casefile = load_caseopfmat(self.cases[i])
            # read case file and return a case dictionary
            case = read_case_from_mat_file(casefile)
            # solve opf model and return a solution and the final objective value
            solution, objval = solve_opf_model(settings, case)
            # check whether the solution point looks correct
            self.assertTrue(solution is not None)
            self.assertTrue(objval is not None)
            self.assertTrue(len(solution) == self.numvars_dc[i])
            self.assertLess(abs(objval - self.objvals_dc[i]), 1e-4)
            self.assertTrue("lincost" in solution.keys())
            self.assertLess(abs(solution["lincost"] - self.lincost_dc[i]), 1e-4)
            self.assertTrue("theta_9" in solution.keys())
            self.assertLess(abs(solution["theta_9"] - self.theta_9_dc[i]), 1e-4)

    # test DC formulation
    def test_dcopf(self):
        # load path to settings file
        settingsfile = load_dcopfsettings()
        # read settings file and return a settings dictionary
        settings = read_settings_from_file(settingsfile)

        for i in range(self.numcases):
            # load path to case file
            casefile = load_caseopf(self.cases[i])
            # read case file and return a case dictionary
            case = read_case_from_file(casefile)
            # solve opf model and return a solution and the final objective value
            solution, objval = solve_opf_model(settings, case)
            # check whether the solution point looks correct
            self.assertTrue(solution is not None)
            self.assertTrue(objval is not None)
            self.assertTrue(len(solution) == self.numvars_dc[i])
            self.assertLess(abs(objval - self.objvals_dc[i]), 1e-4)
            self.assertTrue("lincost" in solution.keys())
            self.assertLess(abs(solution["lincost"] - self.lincost_dc[i]), 1e-4)
            self.assertTrue("theta_9" in solution.keys())
            self.assertLess(abs(solution["theta_9"] - self.theta_9_dc[i]), 1e-4)

    # test AC formulation
    def test_acopf(self):
        # load path to settings file
        settingsfile = load_acopfsettings()
        # read settings file and return a settings dictionary
        settings = read_settings_from_file(settingsfile)
        # load path to case file
        casefile = load_caseopf("9")
        # read case file and return a case dictionary
        case = read_case_from_file(casefile)
        # solve opf model and return a solution and the final objective value
        solution, objval = solve_opf_model(settings, case)
        # check whether the solution points looks correct
        self.assertTrue(solution is not None)
        self.assertTrue(objval is not None)
        self.assertTrue(len(solution) == 134)
        self.assertLess(abs(objval - 5296.665647261), 1e-4)
        self.assertTrue("lincost" in solution.keys())
        self.assertLess(abs(solution["lincost"] - 704.399018), 1e-4)
        self.assertTrue("twinP_1_1_4" in solution.keys())
        self.assertLess(abs(solution["twinP_1_1_4"]), 1e-4)
        self.assertTrue("Q_9_4_9" in solution.keys())
        self.assertLess(abs(solution["Q_9_4_9"] - 0.127732), 1e-4)

    # test IV formulation
    def test_ivopf(self):
        # load path to settings file
        settingsfile = load_ivopfsettings()
        # read settings file and return a settings dictionary
        settings = read_settings_from_file(settingsfile)
        # load path to case file
        casefile = load_caseopf("9")
        # read case file and return a case dictionary
        case = read_case_from_file(casefile)
        # solve opf model and return a solution and the final objective value
        solution, objval = solve_opf_model(settings, case)
        # check whether the solution points looks correct
        self.assertTrue(solution is not None)
        self.assertTrue(objval is not None)
        self.assertTrue(len(solution) == 62)
        self.assertLess(abs(objval - 5296.716652), 1e-4)
        self.assertTrue("lincost" in solution.keys())
        self.assertLess(abs(solution["lincost"] - 703.019340), 1e-4)
        self.assertTrue("f_9" in solution.keys())
        self.assertLess(abs(solution["f_9"] + 0.085954), 1e-4)
        self.assertTrue("P_9_4_9" in solution.keys())
        self.assertLess(abs(solution["P_9_4_9"] - 0.543771), 1e-4)

    def test_graphics(self):
        # load settings dictionary
        settings_graphics_dict = load_opfdictgraphicssettings()
        # load case dictionary
        case = load_opfdictcase()
        # load a precomputed solution and objective value
        solution, objval = load_case9solution()
        # plot the given solution
        plot_opf_solution(settings_graphics_dict, case, solution, objval)

    def test_graphics_after_solving(self):
        # load settings dictionary
        settings_optimization_dict = load_opfdictsettings()
        # load case dictionary
        case = load_opfdictcase()
        # solve opf model and return a solution and the final objective value
        solution, objval = solve_opf_model(settings_optimization_dict, case)
        # plot the computed solution
        settings_graphics_dict = load_opfdictgraphicssettings()
        plot_opf_solution(settings_graphics_dict, case, solution, objval)

    def test_graphics_settings_file(self):
        # load path to settings file
        settingsfile = load_opfgraphicssettings()
        # read settings file and return a settings dictionary
        settings_graphics_dict = read_settings_from_file(settingsfile, True)
        # load path to case file
        casefile = load_caseopf("9")
        # read case file and return a case dictionary
        case = read_case_from_file(casefile)
        # load a precomputed solution and objective value
        solution, objval = load_case9solution()
        # plot the computed solution
        plot_opf_solution(settings_graphics_dict, case, solution, objval)
