import unittest

from gurobi_optimods.opf import solve_acopf_model
from gurobi_optimods.datasets import load_opf


class TestOpf(unittest.TestCase):
    def test_simple(self):
        conf, case = load_opf()
        solution = solve_acopf_model(conf, case)
        self.assertTrue(solution is None)

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
