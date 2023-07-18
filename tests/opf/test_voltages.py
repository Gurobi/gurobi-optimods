# Tests of the voltage/violations public API functions

import unittest

from gurobi_optimods.datasets import load_filepath, load_opfdictcase
from gurobi_optimods.opf import (
    compute_violations_from_given_voltages,
    read_voltages_from_csv_file,
)


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
