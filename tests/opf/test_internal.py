# Tests that solve models using the internal API (more available parameters)

import unittest
from dataclasses import dataclass
from typing import List

import gurobipy as gp
from gurobipy import GRB

from gurobi_optimods.opf import grbformulator_common


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
