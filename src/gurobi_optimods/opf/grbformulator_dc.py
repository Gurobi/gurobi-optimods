import logging
import math

import gurobipy as gp
from gurobipy import GRB

from gurobi_optimods.opf.grbformulator_ac import computebalbounds
from gurobi_optimods.opf.grbformulator_common import set_gencost_objective

logger = logging.getLogger(__name__)


def lpformulator_dc_body(alldata, model):
    """Add variables and constraints for DC formulation to the given model"""

    # Assert compatible settings: Avoiding these asserts when called from the
    # public API should be handled correctly by build_internal_settings(...).
    # 1. branchswitching = 2 only available for AC. Turning it off.
    assert not alldata["branchswitching_comp"]

    _add_dc_gen_bus_variables(alldata, model)
    _add_dc_branch_variables(alldata, model)
    model.update()  # needed to access the variable bounds
    set_gencost_objective(alldata, model)
    _add_dc_branch_activepower_constraints(alldata, model)
    _add_dc_bus_balance_constraints(alldata, model)
    _add_dc_bus_injection_constraints(alldata, model)
    _add_min_active_branch_constraint(alldata, model)


def _add_dc_gen_bus_variables(alldata, model):
    """
    Add variables associated with buses and generators

    - thetavar: voltage angle (voltage magnitude is always 1 for DC)
    - Pinjvar: variable modeling total active power injected by bus j into the
      branches incident with j
    - GenPVar: DC generator real power injection variables
    """

    fixtolerance = alldata["fixtolerance"]
    buses = alldata["buses"]
    gens = alldata["gens"]

    thetavar = {}
    Pinjvar = {}
    GenPvar = {}

    for j, bus in buses.items():
        ubound = 2 * math.pi
        lbound = -ubound
        if bus.inputvoltage:
            lbound = max(lbound, bus.inputA_rad - fixtolerance)
            ubound = min(ubound, bus.inputA_rad + fixtolerance)
        thetavar[bus] = model.addVar(lb=lbound, ub=ubound, name=f"theta_{bus.nodeID}")

        Pubound, Plbound, _, _ = computebalbounds(alldata, bus)
        Pinjvar[bus] = model.addVar(lb=Plbound, ub=Pubound, name=f"IP_{j}")

        # This is done in the inner loop because each generator should have only
        # one associated bus. May be cleaner to just iterate over the generators.
        for genid in bus.genidsbycount:
            gen = gens[genid]
            lower = gen.Pmin * gen.status
            upper = gen.Pmax * gen.status
            # if bus.nodetype == 3:
            #     upper = GRB.INFINITY
            #     lower = -GRB.INFINITY  #ignoring slack bus
            GenPvar[gen] = model.addVar(
                lb=lower, ub=upper, name=f"GP_{gen.count}_{gen.nodeID}"
            )

    alldata["LP"]["thetavar"] = thetavar
    alldata["LP"]["Pinjvar"] = Pinjvar
    alldata["LP"]["GenPvar"] = GenPvar


def _add_dc_branch_variables(alldata, model):
    """
    Add variables associated with branches

    Always:
        - Pvar_f: DC branch real power injected into "from" end of branch. DC
          branch real power injected into "to" end of branch is the same as
          Pvar_f.

    If branch switching is used:
        - twinPvar_f: Auxiliary variable in case branch-switching is being used
        - zvar: branch switching decision variable
    """

    branches = alldata["branches"]
    branchswitching = alldata["branchswitching_mip"]

    Pvar_f = {}  # DC, so f-flow = - t-flow
    twinPvar_f = {}
    zvar = {}

    for j, branch in branches.items():
        bound = branch.limit if branch.constrainedflow else alldata["sumPd"]
        Pvar_f[branch] = model.addVar(
            lb=-bound,
            ub=bound,
            name=f"P_{j}_{branch.f}_{branch.t}",
        )

        if branchswitching:
            twinPvar_f[branch] = model.addVar(
                lb=-GRB.INFINITY,
                ub=GRB.INFINITY,
                name=f"twinP_{j}_{branch.f}_{branch.t}",
            )

    if branchswitching:
        for j, branch in branches.items():
            zvar[branch] = model.addVar(
                obj=0.0,
                vtype=GRB.BINARY,
                name=f"z_{j}_{branch.f}_{branch.t}",
            )

    alldata["LP"]["Pvar_f"] = Pvar_f
    alldata["LP"]["twinPvar_f"] = twinPvar_f
    alldata["MIP"]["zvar"] = zvar


def _add_dc_branch_activepower_constraints(alldata, model):
    """Add constraints defining active power"""

    buses = alldata["buses"]
    branches = alldata["branches"]
    IDtoCountmap = alldata["IDtoCountmap"]
    thetavar = alldata["LP"]["thetavar"]
    Pvar_f = alldata["LP"]["Pvar_f"]
    twinPvar_f = alldata["LP"]["twinPvar_f"]
    zvar = alldata["MIP"]["zvar"]
    branchswitching = alldata["branchswitching_mip"]

    for j, branch in branches.items():
        busf = buses[IDtoCountmap[branch.f]]
        bust = buses[IDtoCountmap[branch.t]]

        if not branch.status:  # out of operation
            branch.Pdeffconstr = model.addConstr(
                Pvar_f[branch] == 0, name=f"Pdef_{j}_{branch.f}_{branch.t}"
            )
            continue

        coeff = 1 / (branch.x * branch.ratio)
        if branchswitching:
            lhs = Pvar_f[branch] + twinPvar_f[branch]
        else:
            lhs = Pvar_f[branch]

        branch.Pdeffconstr = model.addConstr(
            lhs
            == coeff * thetavar[busf]
            - coeff * thetavar[bust]
            - coeff * branch.angle_rad,
            name=f"Pdef_{j}_{branch.f}_{branch.t}",
        )

        if branchswitching:
            coeff = branch.limit if branch.constrainedflow else alldata["sumPd"]

            model.addConstr(
                Pvar_f[branch] <= coeff * zvar[branch],
                name=f"upmip_{j}_{branch.f}_{branch.t}",
            )
            model.addConstr(
                Pvar_f[branch] >= -coeff * zvar[branch],
                name=f"dnmip_{j}_{branch.f}_{branch.t}",
            )

            coeff = 1 / (branch.x * branch.ratio)
            bigM = coeff * max(
                thetavar[busf].UB - thetavar[bust].LB,
                thetavar[bust].UB - thetavar[busf].LB,
            )
            model.addConstr(
                twinPvar_f[branch] <= bigM * (1 - zvar[branch]),
                name=f"upmip_twin_{j}_{branch.f}_{branch.t}",
            )
            model.addConstr(
                twinPvar_f[branch] >= -bigM * (1 - zvar[branch]),
                name=f"dnmip_twin__{j}_{branch.f}_{branch.t}",
            )


def _add_dc_bus_balance_constraints(alldata, model):
    """Balance definintions"""

    buses = alldata["buses"]
    branches = alldata["branches"]
    Pvar_f = alldata["LP"]["Pvar_f"]
    Pinjvar = alldata["LP"]["Pinjvar"]

    alldata["LP"]["balancecons"] = {
        bus: model.addConstr(
            gp.quicksum(
                Pvar_f[branches[branchid]] for branchid in bus.frombranchids.values()
            )
            - gp.quicksum(
                Pvar_f[branches[branchid]] for branchid in bus.tobranchids.values()
            )
            == Pinjvar[bus],
            name=f"PBaldef_{j}_{bus.nodeID}",
        )
        for j, bus in buses.items()
    }


def _add_dc_bus_injection_constraints(alldata, model):
    """Injection definitions"""

    buses = alldata["buses"]
    gens = alldata["gens"]
    Pinjvar = alldata["LP"]["Pinjvar"]
    GenPvar = alldata["LP"]["GenPvar"]

    for j, bus in buses.items():
        model.addConstr(
            Pinjvar[bus]
            == gp.quicksum(GenPvar[gens[genid]] for genid in bus.genidsbycount)
            - bus.Pd,
            name=f"Bus_PInj_{j}",
        )


def _add_min_active_branch_constraint(alldata, model):
    branches = alldata["branches"]
    zvar = alldata["MIP"]["zvar"]
    N = math.floor(alldata["numbranches"] * alldata["minactivebranches"])

    if alldata["branchswitching_mip"]:
        logger.info(f"In bound_zs constraint, {N=}")
        model.addConstr(
            gp.quicksum(zvar[branch] for _, branch in branches.items()) >= N,
            name="sumz_lower_heuristic_bound",
        )
