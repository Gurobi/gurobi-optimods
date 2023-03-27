import unittest

from gurobi_optimods.opf import solve_opf_model
from gurobi_optimods.datasets import load_acopf, load_dcopf


class TestOpf(unittest.TestCase):
    # test simple is on purpose the same as test_acopf for now
    def test_simple(self):
        conf, case = load_acopf()
        solution, objval = solve_opf_model(conf, case)
        # check whether the solution points looks correct
        self.assertTrue(solution is not None)
        self.assertTrue(objval is not None)
        self.assertTrue(len(solution) == 134)
        self.assertLess(abs(objval - 5296.665647261), 1e-4)
        self.assertTrue("lincost" in solution.keys())
        self.assertLess(abs(solution["lincost"] - 704.3990181513464), 1e-4)
        self.assertTrue("twinP_1_1_4" in solution.keys())
        self.assertLess(abs(solution["twinP_1_1_4"]), 1e-9)
        self.assertTrue("Q_9_4_9" in solution.keys())
        self.assertLess(abs(solution["Q_9_4_9"] - 0.12773254991652677), 1e-4)

    def test_acopf(self):
        conf, case = load_acopf()
        solution, objval = solve_opf_model(conf, case)
        # check whether the solution points looks correct
        self.assertTrue(solution is not None)
        self.assertTrue(objval is not None)
        self.assertTrue(len(solution) == 134)
        self.assertLess(abs(objval - 5296.665647261), 1e-4)
        self.assertTrue("lincost" in solution.keys())
        self.assertLess(abs(solution["lincost"] - 704.3990181513464), 1e-4)
        self.assertTrue("twinP_1_1_4" in solution.keys())
        self.assertLess(abs(solution["twinP_1_1_4"]), 1e-9)
        self.assertTrue("Q_9_4_9" in solution.keys())
        self.assertLess(abs(solution["Q_9_4_9"] - 0.12773254991652677), 1e-4)

    def test_dcopf(self):
        conf, case = load_dcopf()
        solution, objval = solve_opf_model(conf, case)
        # check whether the solution point looks correct
        self.assertTrue(solution is not None)
        self.assertTrue(objval is not None)
        self.assertTrue(len(solution) == 50)
        self.assertLess(abs(objval - 5216.026607747), 1e-4)
        self.assertTrue("lincost" in solution.keys())
        self.assertLess(abs(solution["lincost"] - 688.1335088379091), 1e-4)
        self.assertTrue("theta_9" in solution.keys())
        self.assertLess(abs(solution["theta_9"] - 6.083000384352472), 1e-4)
        self.assertTrue("z_9_9_4" in solution.keys())
        self.assertLess(abs(solution["z_9_9_4"] - 1), 1e-4)

    # TODO:
    #  - solve_acopf_model should return something to solution and check it
    #  - add test for plotting
    #  - use logging package from Python instead of own logger
    #  - we should also test dc and iv formulation and solving
    #  - rename files, currently they are called grb<something> or just <something>
    #  - get rid of break_exits
    #  - how to deal with log and other files?
    #    - we probably want to turn off all logs from Gurobi and run everything through a message callback
    #      to avoid a gurobi.log and an OPF.log
    #    - we probably don't want to write any file by default
    #  - there are additional package requirements for plotting
    #  - how do we want to do error handling
    #  - there is a global variable in grbformulator_dc which does not have to be one (it is used in callback)
    #
