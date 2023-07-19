# Tests of plotting functions

import unittest

import gurobipy as gp

from gurobi_optimods.datasets import load_caseNYopf, load_filepath
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


@unittest.skipIf(plotly is None, "plotly is not installed")
class TestGraphicsCase9(unittest.TestCase):
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

        from tests.opf import read_case

        # Info related to case9
        self.case9 = read_case("9")
        self.case9_solution = solve_opf_model(self.case9, opftype="AC", verbose=False)
        self.case9_coords = read_coords_from_csv_file(load_filepath("case9coords.csv"))
        volts_data = read_voltages_from_csv_file(load_filepath("case9volts.csv"))
        self.case9_violations = compute_violations_from_given_voltages(
            self.case9, volts_data, polar=True, verbose=False
        )

        # Manually create a solution with some branches switched off
        from copy import deepcopy

        self.switching_solution = deepcopy(self.case9_solution)
        self.switching_solution["branch"][3]["switching"] = 0
        self.switching_solution["branch"][6]["switching"] = 0

    def test_plot_solution(self):
        # Plot figure using case, coordinates, solution
        fig = generate_opf_solution_figure(
            self.case9, self.case9_coords, self.case9_solution
        )

        # Check whether figure coordinates and scaled input coordinates are the same
        for i in range(9):
            self.assertLess(abs(fig.data[1].x[i] - self.graphics_9_x[i]), 1e-9)
            self.assertLess(abs(fig.data[1].y[i] - self.graphics_9_y[i]), 1e-9)

        # If set to true, plot opens in browser for manual checking
        if False:
            fig.show()

    def test_plot_branchswitching(self):
        # Plot figure using case, coordinates, switching solution
        fig = generate_opf_solution_figure(
            self.case9, self.case9_coords, self.switching_solution
        )

        # If set to true, plot opens in browser for manual checking
        if False:
            fig.show()

    def test_plot_violations(self):
        # Plot violations figure using case, coordinates, voltage solution
        fig = generate_opf_violations_figure(
            self.case9, self.case9_coords, self.case9_violations
        )

        # If set to true, plot opens in browser for manual checking
        if False:
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

        if self.plot_graphics:
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
