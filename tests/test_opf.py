import unittest

from gurobi_optimods.opf import (
    solve_opf_model,
    plot_opf_solution,
    read_settings_from_file,
    read_case_from_file,
)
from gurobi_optimods.datasets import (
    load_case9opf,
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
    # test simple is on purpose the same as test_acopf for now
    # will be removed in final version
    def test_simple(self):
        # load path to settings file
        settingsfile = load_simpleopfsettings()
        # read settings file and return a settings dictionary
        settings = read_settings_from_file(settingsfile)
        # load path to case file
        casefile = load_case9opf()
        # read case file and return a case dictionary
        case = read_case_from_file(casefile)
        # solve opf model and return a solution and the final objective value
        solution, objval = solve_opf_model(settings, case)
        # check whether the solution points looks correct
        self.assertTrue(solution is not None)
        self.assertTrue(objval is not None)

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
        casefile = load_case9opf()
        # read case file and return a case dictionary
        case = read_case_from_file(casefile)
        # load a precomputed solution and objective value
        solution, objval = load_case9solution()
        # plot the computed solution
        plot_opf_solution(settings_graphics_dict, case, solution, objval)

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

    # test AC formulation
    def test_acopf(self):
        # load path to settings file
        settingsfile = load_acopfsettings()
        # read settings file and return a settings dictionary
        settings = read_settings_from_file(settingsfile)
        # load path to case file
        casefile = load_case9opf()
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

    # test DC formulation
    def test_dcopf(self):
        # load path to settings file
        settingsfile = load_dcopfsettings()
        # read settings file and return a settings dictionary
        settings = read_settings_from_file(settingsfile)
        # load path to case file
        casefile = load_case9opf()
        # read case file and return a case dictionary
        case = read_case_from_file(casefile)
        # solve opf model and return a solution and the final objective value
        solution, objval = solve_opf_model(settings, case)
        # check whether the solution point looks correct
        self.assertTrue(solution is not None)
        self.assertTrue(objval is not None)
        self.assertTrue(len(solution) == 50)
        self.assertLess(abs(objval - 5216.026607), 1e-4)
        self.assertTrue("lincost" in solution.keys())
        self.assertLess(abs(solution["lincost"] - 688.133508), 1e-4)
        self.assertTrue("theta_9" in solution.keys())
        self.assertLess(abs(solution["theta_9"] - 6.083), 1e-4)
        self.assertTrue("z_9_9_4" in solution.keys())
        self.assertLess(abs(solution["z_9_9_4"] - 1), 1e-4)

    # test IV formulation
    def test_ivopf(self):
        # load path to settings file
        settingsfile = load_ivopfsettings()
        # read settings file and return a settings dictionary
        settings = read_settings_from_file(settingsfile)
        # load path to case file
        casefile = load_case9opf()
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
