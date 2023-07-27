# Tests that solve models using the internal API (more available parameters)

import os
import unittest
from dataclasses import dataclass
from typing import List

import gurobipy as gp
from gurobipy import GRB

from gurobi_optimods.datasets import load_opf_example
from gurobi_optimods.opf import converters, grbformulator, grbformulator_common
from gurobi_optimods.opf.api import _solve_opf_model_internal


@unittest.skipIf(
    os.environ.get("CI", "false") == "true",
    "Tests expensive internal options; skipped for CI",
)
class TestInternal(unittest.TestCase):
    # Test internal options we haven't exposed on the public API yet

    def setUp(self):
        self.env = gp.Env()
        self.case = load_opf_example("case9")

    def tearDown(self):
        self.env.close()

    def test_ac_polar(self):
        # Test polar formulation. It's expensive to test all combinations,
        # so for now just test the defaults.
        self.env.setParam("SolutionLimit", 1)
        solution = _solve_opf_model_internal(
            self.env,
            self.case,
            opftype="AC",
            polar=True,
            useef=True,
            usejabr=True,
            ivtype="aggressive",
            branchswitching=0,
            usemipstart=True,
            minactivebranches=0.9,
            useactivelossineqs=False,
        )
        self.assertIsNotNone(solution)
        self.assertEqual(solution["success"], 1)

    def test_ac_settings(self):
        # Test all AC solver options with polar=False.
        settingslist = [
            dict(
                opftype="AC",
                polar=False,
                useef=useef,
                usejabr=usejabr,
                branchswitching=branchswitching,
                usemipstart=usemipstart,
                useactivelossineqs=useactivelossineqs,
                minactivebranches=minactivebranches,
                ivtype="aggressive",
            )
            for useef in [False, True]
            for usejabr in [False, True]
            for branchswitching in [0, 1, 2]
            for usemipstart in [False, True]
            for useactivelossineqs in [False, True]
            for minactivebranches in [0, 0.5, 0.95]
        ]

        # Solve opf model for the same case and return a solution.
        # The model has to be feasible for every setting combination.
        self.env.setParam("SolutionLimit", 1)
        for acopf_settings in settingslist:
            with self.subTest(acopf_settings):
                solution = _solve_opf_model_internal(
                    self.env,
                    self.case,
                    **acopf_settings,
                )
                self.assertIsNotNone(solution)
                self.assertEqual(solution["success"], 1)

    def test_dc_settings(self):
        # DC has more limited configurations: polar, ef, jabr,
        # branchswitching=2, and ivtype are not available
        settingslist = [
            dict(
                opftype="DC",
                branchswitching=branchswitching,
                usemipstart=usemipstart,
                useactivelossineqs=useactivelossineqs,
                minactivebranches=minactivebranches,
                polar=False,
                useef=True,
                usejabr=True,
                ivtype="aggressive",
            )
            for branchswitching in [0, 1]
            for usemipstart in [False, True]
            for useactivelossineqs in [False, True]
            for minactivebranches in [0, 0.5, 0.95]
        ]

        # Solve opf model for the same case and return a solution.
        # The model has to be feasible for every setting combination.
        self.env.setParam("SolutionLimit", 1)
        for acopf_settings in settingslist:
            with self.subTest(acopf_settings):
                solution = _solve_opf_model_internal(
                    self.env,
                    self.case,
                    **acopf_settings,
                )
                self.assertIsNotNone(solution)
                self.assertEqual(solution["success"], 1)

    def test_iv_settings(self):
        # IV has more limited configurations: polar, ef, jabr,
        # branchswitching, and mipstart are not available
        settingslist = [
            dict(
                opftype="IV",
                ivtype=ivtype,
                useactivelossineqs=useactivelossineqs,
                minactivebranches=minactivebranches,
                polar=False,
                useef=True,
                usejabr=True,
                branchswitching=0,
                usemipstart=True,
            )
            for ivtype in ["plain", "aggressive"]
            for useactivelossineqs in [False, True]
            for minactivebranches in [0, 0.5, 0.95]
        ]

        # Solve opf model for the same case and return a solution.
        # The model has to be feasible for every setting combination.
        self.env.setParam("SolutionLimit", 1)
        for acopf_settings in settingslist:
            with self.subTest(acopf_settings):
                solution = _solve_opf_model_internal(
                    self.env,
                    self.case,
                    **acopf_settings,
                )
                self.assertIsNotNone(solution)
                self.assertEqual(solution["success"], 1)

    # test IV formulation
    def test_ivopf(self):
        # currently all other cases take very long in IV formulation
        # solve opf model and return a solution
        solution = _solve_opf_model_internal(
            self.env,
            self.case,
            opftype="IV",
            ivtype="aggressive",
            polar=False,
            useef=True,
            usejabr=True,
            branchswitching=0,
            usemipstart=True,
            minactivebranches=0.9,
            useactivelossineqs=False,
        )
        # check whether the solution points looks correct
        self.assertIsNotNone(solution)
        self.assertEqual(solution["success"], 1)
        self.assertIsNotNone(solution["f"])
        # differences can be quite big because we solve only to 0.1% optimality
        self.assertLess(abs(solution["f"] - 5297.0142014), 1e1)
        self.assertLess(abs(solution["bus"][0]["Vm"] - 1.0986917), 1e-1)
        self.assertLess(abs(solution["gen"][1]["Pg"] - 132.87813), 1e1)
        self.assertLess(abs(solution["gen"][2]["Qg"] + 22.802347), 1e1)
        self.assertLess(abs(solution["branch"][3]["Pf"] - 95.113306), 1e1)
        self.assertLess(abs(solution["branch"][4]["Qt"] + 18.431373), 1e1)


class TestFingerprints(unittest.TestCase):
    # Check model fingerprints for some specific cases. Useful while
    # refactoring the code, but

    def setUp(self):
        self.env = gp.Env()
        self.data = [
            {
                "name": "case9",
                "case": load_opf_example("case9"),
                "dc_fingerprint": 1980532444,
                "dc_switching_fingerprint": 1441352667,
                "ac_fingerprint": -121125607,
                "ac_relax_fingerprint": -552798165,
                "iv_fingerprint": 696846664,
            },
            {
                "name": "caseNY",
                "case": load_opf_example("caseNY"),
                "dc_fingerprint": 1466917423,
                "dc_switching_fingerprint": 1343022114,
                "ac_fingerprint": 1587606503,
                "ac_relax_fingerprint": -557544349,
                "iv_fingerprint": -199906229,
            },
        ]

    def tearDown(self):
        self.env.close()

    def test_dc(self):
        for example in self.data:
            with self.subTest(name=example["name"]):
                with gp.Model(env=self.env) as model:
                    alldata = converters.convert_case_to_internal_format(
                        example["case"]
                    )
                    alldata.update(
                        converters.build_internal_settings(
                            opftype="DC",
                            branchswitching=False,
                            usemipstart=False,
                            useactivelossineqs=False,
                            minactivebranches=0.0,
                            polar=False,
                            useef=True,
                            usejabr=True,
                            ivtype="aggressive",
                        )
                    )
                    grbformulator.lpformulator_setup(alldata, grbformulator.OpfType.DC)
                    grbformulator.lpformulator_body(
                        alldata, model, grbformulator.OpfType.DC
                    )
                    model.update()
                    self.assertEqual(model.Fingerprint, example["dc_fingerprint"])

    def test_dc_switching(self):
        for example in self.data:
            with self.subTest(name=example["name"]):
                with gp.Model(env=self.env) as model:
                    alldata = converters.convert_case_to_internal_format(
                        example["case"]
                    )
                    alldata.update(
                        converters.build_internal_settings(
                            opftype="DC",
                            branchswitching=True,
                            usemipstart=False,
                            useactivelossineqs=False,
                            minactivebranches=0.9,
                            polar=False,
                            useef=True,
                            usejabr=True,
                            ivtype="aggressive",
                        )
                    )
                    grbformulator.lpformulator_setup(alldata, grbformulator.OpfType.DC)
                    grbformulator.lpformulator_body(
                        alldata, model, grbformulator.OpfType.DC
                    )
                    model.update()
                    self.assertEqual(
                        model.Fingerprint, example["dc_switching_fingerprint"]
                    )

    def test_ac(self):
        for example in self.data:
            with self.subTest(name=example["name"]):
                with gp.Model(env=self.env) as model:
                    alldata = converters.convert_case_to_internal_format(
                        example["case"]
                    )
                    alldata.update(
                        converters.build_internal_settings(
                            opftype="AC",
                            polar=False,
                            useef=True,
                            usejabr=True,
                            branchswitching=False,
                            usemipstart=False,
                            useactivelossineqs=False,
                            minactivebranches=0.0,
                            ivtype="aggressive",
                        )
                    )
                    grbformulator.lpformulator_setup(alldata, grbformulator.OpfType.AC)
                    grbformulator.lpformulator_body(
                        alldata, model, grbformulator.OpfType.AC
                    )
                    model.update()
                    self.assertEqual(model.Fingerprint, example["ac_fingerprint"])

    def test_ac_relax(self):
        for example in self.data:
            with self.subTest(name=example["name"]):
                with gp.Model(env=self.env) as model:
                    alldata = converters.convert_case_to_internal_format(
                        example["case"]
                    )
                    alldata.update(
                        converters.build_internal_settings(
                            opftype="AC",
                            polar=False,
                            useef=False,
                            usejabr=True,
                            branchswitching=False,
                            usemipstart=False,
                            useactivelossineqs=False,
                            minactivebranches=0.0,
                            ivtype="aggressive",
                        )
                    )
                    grbformulator.lpformulator_setup(alldata, grbformulator.OpfType.AC)
                    grbformulator.lpformulator_body(
                        alldata, model, grbformulator.OpfType.AC
                    )
                    model.update()
                    self.assertEqual(model.Fingerprint, example["ac_relax_fingerprint"])

    def test_iv(self):
        for example in self.data:
            with self.subTest(name=example["name"]):
                with gp.Model(env=self.env) as model:
                    alldata = converters.convert_case_to_internal_format(
                        example["case"]
                    )
                    alldata.update(
                        converters.build_internal_settings(
                            opftype="IV",
                            ivtype="aggressive",
                            useactivelossineqs=False,
                            minactivebranches=0.0,
                            polar=False,
                            useef=True,
                            usejabr=True,
                            branchswitching=0,
                            usemipstart=True,
                        )
                    )
                    grbformulator.lpformulator_setup(alldata, grbformulator.OpfType.IV)
                    grbformulator.lpformulator_body(
                        alldata, model, grbformulator.OpfType.IV
                    )
                    model.update()
                    self.assertEqual(model.Fingerprint, example["iv_fingerprint"])


@dataclass
class Gen:
    costdegree: int
    costvector: List[float]

    def __hash__(self):
        return id(self)


class TestGenCostObjective(unittest.TestCase):
    # Small test cases of the objective formulation. Compare against reference
    # formulation.

    def setUp(self):
        self.env = gp.Env()

    def tearDown(self):
        self.env.close()

    def test_linear(self):
        # Two generators with linear cost

        gen1 = Gen(costdegree=1, costvector=[2.0, 3.0])
        gen2 = Gen(costdegree=1, costvector=[5.0, 6.0])

        with gp.Model(env=self.env) as refmodel:
            x1, x2 = refmodel.addVars([1, 2], name="x").values()
            lincostvar = refmodel.addVar(lb=-GRB.INFINITY, name="lincost")
            constvar = refmodel.addVar(lb=1.0, ub=1.0, name="constant")
            refmodel.addConstr(2.0 * x1 + 5.0 * x2 == lincostvar, name="lincostdef")
            refmodel.setObjective(9.0 * constvar + lincostvar)
            refmodel.update()
            ref_fingerprint = refmodel.Fingerprint

        with gp.Model(env=self.env) as model:
            x1, x2 = model.addVars([1, 2], name="x").values()
            alldata = {
                "LP": {"GenPvar": {gen1: x1, gen2: x2}},
                "gens": {1: gen1, 2: gen2},
            }

            grbformulator_common.set_gencost_objective(alldata, model)
            model.update()
            self.assertEqual(model.Fingerprint, ref_fingerprint)

    def test_quadratic(self):
        # Two generators with quadratic cost

        gen1 = Gen(costdegree=2, costvector=[2.0, 3.0, 4.0])
        gen2 = Gen(costdegree=2, costvector=[5.0, 6.0, 7.0])

        with gp.Model(env=self.env) as refmodel:
            x1, x2 = refmodel.addVars([1, 2], name="x").values()
            lincostvar = refmodel.addVar(lb=-GRB.INFINITY, name="lincost")
            constvar = refmodel.addVar(lb=1.0, ub=1.0, name="constant")
            refmodel.addConstr(3.0 * x1 + 6.0 * x2 == lincostvar, name="lincostdef")
            refmodel.setObjective(
                11.0 * constvar + lincostvar + 2.0 * x1 * x1 + 5.0 * x2 * x2
            )
            refmodel.update()
            ref_fingerprint = refmodel.Fingerprint

        with gp.Model(env=self.env) as model:
            x1, x2 = model.addVars([1, 2], name="x").values()
            alldata = {
                "LP": {"GenPvar": {gen1: x1, gen2: x2}},
                "gens": {1: gen1, 2: gen2},
            }

            grbformulator_common.set_gencost_objective(alldata, model)
            model.update()
            self.assertEqual(model.Fingerprint, ref_fingerprint)

    def test_degree3_quadratic(self):
        # Two generators with quadratic cost, but costdegree 3 (cubic terms zero)

        gen1 = Gen(costdegree=3, costvector=[0.0, 2.0, 3.0, 4.0])
        gen2 = Gen(costdegree=3, costvector=[0.0, 5.0, 6.0, 7.0])

        with gp.Model(env=self.env) as refmodel:
            x1, x2 = refmodel.addVars([1, 2], name="x").values()
            lincostvar = refmodel.addVar(lb=-GRB.INFINITY, name="lincost")
            constvar = refmodel.addVar(lb=1.0, ub=1.0, name="constant")
            refmodel.addConstr(3.0 * x1 + 6.0 * x2 == lincostvar, name="lincostdef")
            refmodel.setObjective(
                11.0 * constvar + lincostvar + 2.0 * x1 * x1 + 5.0 * x2 * x2
            )
            refmodel.update()
            ref_fingerprint = refmodel.Fingerprint

        with gp.Model(env=self.env) as model:
            x1, x2 = model.addVars([1, 2], name="x").values()
            alldata = {
                "LP": {"GenPvar": {gen1: x1, gen2: x2}},
                "gens": {1: gen1, 2: gen2},
            }

            grbformulator_common.set_gencost_objective(alldata, model)
            model.update()
            self.assertEqual(model.Fingerprint, ref_fingerprint)
