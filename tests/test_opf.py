import unittest

from gurobi_optimods.opf import solve_opf_model, plot_opf_solution
from gurobi_optimods.datasets import (
    load_case9opf,
    load_acopfsettings,
    load_dcopfsettings,
    load_ivopfsettings,
    load_simpleopfsettings,
    load_opfdictcase,
    load_opfdictsettings,
    load_opfdictgraphicssettings,
    load_case9solution,
)


class TestOpf(unittest.TestCase):
    # test simple is on purpose the same as test_acopf for now
    # will be removed in final version
    def test_simple(self):
        conf = load_simpleopfsettings()
        case = load_case9opf()
        solution, objval = solve_opf_model(conf, case)
        # check whether the solution points looks correct
        self.assertTrue(solution is not None)
        self.assertTrue(objval is not None)

    def test_graphics(self):
        conf_graphics = load_opfdictgraphicssettings()
        case = load_opfdictcase()
        solution, objval = load_case9solution()
        plot_opf_solution(conf_graphics, case, solution, objval)

    def test_graphics_after_solving(self):
        conf_optimization = load_opfdictsettings()
        case = load_opfdictcase()
        solution, objval = solve_opf_model(conf_optimization, case)
        conf_graphics = load_opfdictgraphicssettings()
        plot_opf_solution(conf_graphics, case, solution, objval)

    # test reading case from a dict
    def test_opfdictcase(self):
        conf = load_acopfsettings()
        case = load_opfdictcase()
        solution, objval = solve_opf_model(conf, case)
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

    # test reading settings from a dict
    def test_opfdictsettings(self):
        conf = load_opfdictsettings()
        case = load_case9opf()
        solution, objval = solve_opf_model(conf, case)
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

    # test reading settings an case file from dicts
    def test_opfdicts(self):
        conf = load_opfdictsettings()
        case = load_opfdictcase()
        solution, objval = solve_opf_model(conf, case)
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
        conf = load_acopfsettings()
        case = load_case9opf()
        solution, objval = solve_opf_model(conf, case)
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
        conf = load_dcopfsettings()
        case = load_case9opf()
        solution, objval = solve_opf_model(conf, case)
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
        conf = load_ivopfsettings()
        case = load_case9opf()
        solution, objval = solve_opf_model(conf, case)
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

    # TODO:
    #  - add test for plotting
    #  - rename files, currently they are called grb<something> and some have own names
    #  - get rid of break_exits
    #  - how to deal with log and other files?
    #    - which files do we want by default
    #  - there are additional package requirements for plotting
    #    - in particular there is this graphviz 3rd party software
    #    - it looks like there is also a pipy graphviz package https://pypi.org/project/graphviz/ Maybe we can use this one?
    #    - we could check whether the graphviz package is installed
    #  - how do we want to do error handling
    #  - there is a global variable in grbformulator_dc which does not have to be one (it is used in callback)
    #
