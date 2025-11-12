# Tests of plotting functions
# FIXME: "unclosed socket" warnings when showing the plot

import gzip
import json
import pathlib
import unittest

from gurobipy import GRB

from gurobi_optimods.datasets import load_opf_example, load_opf_extra
from gurobi_optimods.opf import (
    compute_violations,
    solution_plot,
    solve_opf,
    violation_plot,
)

from ..utils import size_limited_license

# If plotly is not installed, tests will be skipped
try:
    import plotly
except ImportError:
    plotly = None


@unittest.skipIf(GRB.VERSION_MAJOR < 12, "Needs Gurobi 12")
@unittest.skipIf(plotly is None, "plotly is not installed")
class TestGraphicsCase9(unittest.TestCase):
    def setUp(self):
        self.plot_graphics = False
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
        self.case9_solution = solve_opf(self.case9, opftype="AC", verbose=False)
        self.case9_coords = load_opf_extra("case9-coordinates")
        volts_data = load_opf_extra("case9-voltages")
        self.case9_violations = compute_violations(
            self.case9, volts_data, polar=False, verbose=False
        )

        # Load manually created solution with some branches switched off
        self.case9_switching = load_opf_example("case9-switching")
        self.switching_solution = json.loads(
            pathlib.Path(__file__)
            .parent.joinpath("data/case9_switching_solution.json")
            .read_text()
        )

    def test_solution_plot(self):
        # Plot figure using case, coordinates, solution
        fig = solution_plot(self.case9, self.case9_coords, self.case9_solution)

        # Check whether figure coordinates and scaled input coordinates are the same
        for i in range(9):
            self.assertLess(abs(fig.data[1].x[i] - self.graphics_9_x[i]), 1e-9)
            self.assertLess(abs(fig.data[1].y[i] - self.graphics_9_y[i]), 1e-9)

        # If set to true, plot opens in browser for manual checking
        if self.plot_graphics:
            fig.show()

    def test_plot_branchswitching(self):
        # Plot figure using case, coordinates, switching solution
        fig = solution_plot(
            self.case9_switching, self.case9_coords, self.switching_solution
        )

        # If set to true, plot opens in browser for manual checking
        if self.plot_graphics:
            fig.show()

    def test_violation_plot(self):
        # Plot violations figure using case, coordinates, voltage solution
        fig = violation_plot(self.case9, self.case9_coords, self.case9_violations)

        # If set to true, plot opens in browser for manual checking
        if self.plot_graphics:
            fig.show()


@unittest.skipIf(plotly is None, "plotly is not installed")
@unittest.skipIf(size_limited_license(), "size-limited-license")
class TestGraphicsNewYork(unittest.TestCase):
    # test a real data set for New York

    def setUp(self):
        self.plot_graphics = False
        self.case = load_opf_example("caseNY")
        self.coords = load_opf_extra("caseNY-coordinates")
        with gzip.open(
            pathlib.Path(__file__).parent.joinpath(
                "data/ny_dc_switching_solution.json.gz"
            )
        ) as infile:
            self.switching_solution = json.load(infile)

    def test_dc_solution(self):
        # Solve and plot DC solution
        solution = solve_opf(self.case, opftype="DC")
        fig = solution_plot(self.case, self.coords, solution)

        # Test a few coordinates
        self.assertLess(abs(fig.data[1].x[0] - 1381.2), 1e-9)
        self.assertLess(abs(fig.data[1].y[0] - 1203.5), 1e-9)
        self.assertLess(abs(fig.data[1].x[-1] - 837.2), 1e-9)
        self.assertLess(abs(fig.data[1].y[-1] - 511.85), 1e-9)

        # If set to true, plot opens in browser for manual checking
        if self.plot_graphics:
            fig.show()

    def test_branchswitching(self):
        # Plot a pre-loaded DC branch switching solution
        fig = solution_plot(self.case, self.coords, self.switching_solution)

        # If set to true, plot opens in browser for manual checking
        if self.plot_graphics:
            fig.show()
