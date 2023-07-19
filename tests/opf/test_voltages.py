# Tests of the voltage/violations public API functions

import unittest

from gurobi_optimods.opf import compute_violations_from_given_voltages


class TestComputeVoltages(unittest.TestCase):
    # test violation computation out of pre-defined voltage data

    def setUp(self):
        from gurobi_optimods.datasets import load_filepath
        from gurobi_optimods.opf import read_voltages_from_csv_file
        from tests.opf import read_case

        self.case = read_case("9")
        self.volts_data = read_voltages_from_csv_file(load_filepath("case9volts.csv"))

    def test_volts(self):
        violations = compute_violations_from_given_voltages(
            self.case, self.volts_data, polar=True
        )

        self.assertTrue(violations is not None)
        self.assertLess(abs(violations["bus"][1]["Vmviol"]), 1e-9)
        self.assertLess(abs(violations["bus"][3]["Pviol"] + 318.899783), 1e-4)
        self.assertLess(abs(violations["bus"][6]["Qviol"] - 9.8410206), 1e-4)
        self.assertLess(abs(violations["branch"][6]["limitviol"] - 66.33435), 1e-4)
