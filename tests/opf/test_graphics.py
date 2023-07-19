# Tests of plotting functions
# FIXME: "unclosed socket" warnings when showing the plot

import json
import pathlib
import unittest

import gurobipy as gp

from gurobi_optimods.datasets import load_filepath, load_opf_example
from gurobi_optimods.opf import compute_violations_from_given_voltages, solve_opf_model
from gurobi_optimods.opf.io import (
    read_coords_from_csv_file,
    read_voltages_from_csv_file,
)

# If plotly is not installed, some tests will be skipped
try:
    import plotly
except ImportError:
    plotly = None

# If plotly is installed, the opfgraphics module should import ok
if plotly:
    from gurobi_optimods.opf.graphics import plot_solution, plot_violations


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

        # Info related to case9
        self.case9 = load_opf_example("case9")
        self.case9_solution = solve_opf_model(self.case9, opftype="AC", verbose=False)
        self.case9_coords = read_coords_from_csv_file(load_filepath("case9coords.csv"))
        volts_data = read_voltages_from_csv_file(load_filepath("case9volts.csv"))
        self.case9_violations = compute_violations_from_given_voltages(
            self.case9, volts_data, polar=True, verbose=False
        )

        # Load manually created solution with some branches switched off
        self.switching_solution = json.loads(
            pathlib.Path(__file__)
            .parent.joinpath("data/case9_switching_solution.json")
            .read_text()
        )

    def test_plot_solution(self):
        # Plot figure using case, coordinates, solution
        fig = plot_solution(self.case9, self.case9_coords, self.case9_solution)

        # Check whether figure coordinates and scaled input coordinates are the same
        for i in range(9):
            self.assertLess(abs(fig.data[1].x[i] - self.graphics_9_x[i]), 1e-9)
            self.assertLess(abs(fig.data[1].y[i] - self.graphics_9_y[i]), 1e-9)

        # If set to true, plot opens in browser for manual checking
        if False:
            fig.show()

    def test_plot_branchswitching(self):
        # Plot figure using case, coordinates, switching solution
        fig = plot_solution(self.case9, self.case9_coords, self.switching_solution)

        # If set to true, plot opens in browser for manual checking
        if False:
            fig.show()

    def test_plot_violations(self):
        # Plot violations figure using case, coordinates, voltage solution
        fig = plot_violations(self.case9, self.case9_coords, self.case9_violations)

        # If set to true, plot opens in browser for manual checking
        if False:
            fig.show()


@unittest.skipIf(plotly is None, "plotly is not installed")
@unittest.skipIf(size_limited_license(), "size-limited-license")
class TestGraphicsNewYork(unittest.TestCase):
    # test a real data set for New York

    def setUp(self):
        self.case = load_opf_example("caseNY")
        coordsfile = load_filepath("nybuses.csv")
        self.coords = read_coords_from_csv_file(coordsfile)
        self.switching_solution = json.loads(
            pathlib.Path(__file__)
            .parent.joinpath("data/ny_dc_switching_solution.json")
            .read_text()
        )

    def test_dc_solution(self):
        # Solve and plot DC solution
        solution = solve_opf_model(self.case, opftype="DC")
        fig = plot_solution(self.case, self.coords, solution)

        # Test a few coordinates
        self.assertLess(abs(fig.data[1].x[0] - 1381.2), 1e-9)
        self.assertLess(abs(fig.data[1].y[0] - 1203.5), 1e-9)
        self.assertLess(abs(fig.data[1].x[-1] - 837.2), 1e-9)
        self.assertLess(abs(fig.data[1].y[-1] - 511.85), 1e-9)

        # If set to true, plot opens in browser for manual checking
        if False:
            fig.show()

    def test_branchswitching(self):
        # Plot a pre-loaded DC branch switching solution
        fig = plot_solution(self.case, self.coords, self.switching_solution)

        # If set to true, plot opens in browser for manual checking
        if False:
            fig.show()
