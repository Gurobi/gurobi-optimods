import logging
import math
import time
from enum import Enum

import gurobipy as gp
from gurobipy import GRB

from gurobi_optimods.opf.grbformulator_ac import lpformulator_ac_body
from gurobi_optimods.opf.grbformulator_dc import lpformulator_dc_body
from gurobi_optimods.opf.grbformulator_iv import lpformulator_iv_body

logger = logging.getLogger(__name__)


class OpfType(Enum):
    """
    Defines possible OPF formulation types
    """

    AC = "AC"
    DC = "DC"
    IV = "IV"


def construct_and_solve_model(env, alldata):
    """Construct OPF model for given data and solve it. Return the solution
    dictionary."""

    if alldata["doac"]:
        opftype = OpfType.AC
    elif alldata["dodc"]:
        opftype = OpfType.DC
    elif alldata["doiv"]:
        opftype = OpfType.IV
    else:
        # Should not happen (this would be a configuration setup bug)
        raise ValueError("No model type set")

    logger.info(f"Running {opftype.value}OPF formulation.")

    # Create model
    with gp.Model(f"{opftype.value}_Formulation_Model", env=env) as model:
        # Formulate model, update and print statistics

        tic = time.perf_counter()

        if opftype == OpfType.AC:
            lpformulator_ac_body(alldata, model)
        elif opftype == OpfType.DC:
            lpformulator_dc_body(alldata, model)
        elif opftype == OpfType.IV:
            lpformulator_iv_body(alldata, model)
        else:
            raise ValueError("Unknown OPF type.")

        toc = time.perf_counter()

        logger.info(f"{opftype.value}OPF model constructed ({toc - tic:.2f}s). Stats:")
        model.update()
        model.printStats()

        # Solve the OPF model, return a dictionary following MATPOWER notation
        lpformulator_optimize(alldata, model, opftype)
        return turn_solution_into_result_dict(alldata, model, opftype, "result")


def lpformulator_optimize(alldata, model, opftype):
    """
    Optimizes constructed OPF model.

    Resolves model with DualReductions=0 when model is found to
    be infeasible or unbounded in order to get more information.

    Resolves model with additional settings if numerical trouble
    has been encountered to possibly still get a solution.

    Returns the feasible solution count.
    """

    # NonConvex solver needed for normal AC models
    if opftype != OpfType.DC and (alldata["use_ef"] or alldata["dopolar"]):
        model.params.NonConvex = 2

    # Always use a pre-defined MIPStart for DC if we have binary variables
    # For AC only use it if it is requested
    # Note that IV currently does not support branch switching
    if (alldata["branchswitching_mip"] and opftype == OpfType.DC) or (
        alldata["usemipstart"]
        and (alldata["branchswitching_mip"] or alldata["branchswitching_comp"])
    ):
        logger.info("Using mip start with all branches kept on.")
        # MIP start
        # Turn on all branches and hope for the best
        zvar = alldata["MIP"]["zvar"]
        branches = alldata["branches"]
        numbranches = alldata["numbranches"]
        for j in range(1, 1 + numbranches):
            branch = branches[j]
            zvar[branch].Start = 1.0

    model.optimize()

    # Check model status and re-optimize if numerical trouble or inconclusive results.
    if model.status == GRB.INF_OR_UNBD:
        logger.info("Model Status: infeasible or unbounded.")
        logger.info("Re-optimizing with DualReductions turned off.")
        model.Params.DualReductions = 0
        model.optimize()

    if model.status == GRB.INFEASIBLE:
        raise ValueError("Infeasible model")

    if model.status == GRB.NUMERIC:
        logger.info("Solve failed due to numerical issues.")
        logger.info(
            "Re-optimizing with settings focused on improving numerical stability."
        )
        model.Params.NumericFocus = 2
        model.Params.BarHomogeneous = 1
        model.reset()
        model.optimize()

    # Only print objective value and solution quality if at least
    # one feasible point is available
    if model.SolCount > 0:
        logger.info(f"Objective value = {model.objVal}. Solution quality:")
        model.printQuality()

    return model.SolCount


def turn_solution_into_result_dict(alldata, model, opftype, type):
    """
    Turns a Gurobi solution into an OPF dictionary in MATPOWER notation.
    If no solution is present, the result dictionary "success" value is 0.

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param model: Gurobi model which was recently optimized
    :type model: :class: `gurobipy.Model`
    :param opftype: Type of OPF formulation
    :type opftype: :enum: `OpfType`
    :param type: Type of output dictionary.
                 Available are currently "result" and "violations"
    :type type: str

    :raises ValueError: Unknown result type

    :return: Dictionary holding OPF result information following MATPOWER notation
             The "success" entry states whether a feasible solution has been found
             during the optimization process
    :rtype: dict

    TODO most of this could be done by copying the case input dictionary, instead
    of reverse engineering it here
    """

    # Reconstruct case data from our data
    result = {}
    result["baseMVA"] = alldata["baseMVA"]
    baseMVA = result["baseMVA"]
    result["bus"] = {}
    result["gen"] = {}
    result["branch"] = {}
    result["gencost"] = {}

    # Buses
    buses = alldata["buses"]
    index = 1
    for b in buses.values():
        # bus_i type    Pd  Qd  Gs  Bs  area    Vm  Va  baseKV  zone    Vmax    Vmin
        # We don't use area and zone but we save them for consistency with MATPOWER
        matbus = {
            "bus_i": b.nodeID,
            "type": b.nodetype,
            "Pd": b.Pd * baseMVA,
            "Qd": b.Qd * baseMVA,
            "Gs": b.Gs * baseMVA,
            "Bs": b.Bs * baseMVA,
            "area": b.area,
            "Vm": b.Vm,
            "Va": b.Va,
            "baseKV": b.Vbase,
            "zone": b.zone,
            "Vmax": b.Vmax,
            "Vmin": b.Vmin,
        }
        result["bus"][index] = matbus
        index += 1

    # Generators and gen costs
    gens = alldata["gens"]
    index = 1
    for g in gens.values():
        # Generator data
        # bus   Pg  Qg  Qmax    Qmin    Vg  mBase   status  Pmax    Pmin    Pc1 Pc2 Qc1min  Qc1max  Qc2min  Qc2max  ramp_agc    ramp_10 ramp_30 ramp_q  apf
        # We don't use Vg, mBase, and everything after Pmin but we save them for consistency with MATPOWER
        matgen = {
            "bus": g.nodeID,
            "Pg": g.Pg,
            "Qg": g.Qg,
            "Qmax": g.Qmax * baseMVA,
            "Qmin": g.Qmin * baseMVA,
            "Vg": g.Vg,
            "mBase": g.mBase,
            "status": g.status,
            "Pmax": g.Pmax * baseMVA,
            "Pmin": g.Pmin * baseMVA,
            "Pc1": g.Pc1,
            "Pc2": g.Pc2,
            "Qc1min": g.Qc1min,
            "Qc1max": g.Qc1max,
            "Qc2min": g.Qc2min,
            "Qc2max": g.Qc2max,
            "ramp_agc": g.ramp_agc,
            "ramp_10": g.ramp_10,
            "ramp_30": g.ramp_30,
            "ramp_q": g.ramp_q,
            "apf": g.apf,
        }
        result["gen"][index] = matgen
        gencost = g.costvector
        for j in range(len(gencost)):  # scale cost back to MATPOWER format
            gencost[j] /= baseMVA ** (g.costdegree - j)

        # Generator cost data
        #   1   startup shutdown    n   x1  y1  ... xn  yn
        #   2   startup shutdown    n   c(n-1)  ... c0
        # We don't use startup and shutdown but we save them for consistency with MATPOWER
        matgencost = {
            "costtype": g.costtype,
            "startup": g.startup,
            "shutdown": g.shutdown,
            "n": g.costdegree + 1,
            "costvector": gencost,
        }
        result["gencost"][index] = matgencost
        index += 1

    # Branches
    branches = alldata["branches"]
    index = 1
    for b in branches.values():
        # fbus  tbus    r   x   b   rateA   rateB   rateC   ratio   angle   status  angmin  angmax
        matbranch = {
            "fbus": b.f,
            "tbus": b.t,
            "r": b.r,
            "x": b.x,
            "b": b.bc,
            "rateA": b.rateAmva * baseMVA,
            "rateB": b.rateBmva * baseMVA,
            "rateC": b.rateCmva * baseMVA,
            "ratio": b.ratio,
            "angle": b.angle,
            "status": b.status,
            "angmin": b.minangle,
            "angmax": b.maxangle,
        }
        result["branch"][index] = matbranch
        index += 1

    if type == "result":
        fill_result_fields(alldata, model, opftype, result)
    elif type == "violation":
        fill_violations_fields(alldata, opftype, result)
    else:
        raise ValueError("Unknown result type.")

    # Post-conversion to sane python structure
    for field in ["bus", "branch", "gen", "gencost"]:
        entry = result[field]
        result[field] = [entry[i + 1] for i in range(len(entry))]

    return result


def fill_result_fields(alldata, model, opftype, result):
    """
    Extends the result dictionary by additional information derived from Gurobi solution
    See MATPOWER-manual for more details https://matpower.org/docs/MATPOWER-manual.pdf

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param model: Gurobi model which was recently optimized
    :type model: :class: `gurobipy.Model`
    :param opftype: Type of OPF formulation
    :type opftype: :enum: `OpfType`
    :param result: Dictionary holding all case and solution data
    :type result: dict
    """
    baseMVA = result["baseMVA"]

    result["et"] = model.Runtime
    if model.SolCount < 1:
        # We did not find any feasible solution
        result["success"] = 0
        return

    # We have at least one feasible solution
    result["success"] = 1
    result["f"] = model.ObjVal
    # Set DC solution data
    if opftype == OpfType.DC:
        # Bus solution data
        for busindex in result["bus"]:
            resbus = result["bus"][busindex]
            databus = alldata["buses"][busindex]
            # Override old values
            resbus["Va"] = alldata["LP"]["thetavar"][databus].X  # Voltage angle
            resbus["Vm"] = 1  # Voltage magnitude is always 1 for DC
            if not alldata["branchswitching_mip"]:
                resbus["mu"] = alldata["LP"]["balancecons"][
                    databus
                ].Pi  # Shadow prices of balance constraints only available for LP

        # Generator solution data
        for genindex in result["gen"]:
            resgen = result["gen"][genindex]
            datagen = alldata["gens"][genindex]
            # Override old values
            resgen["Pg"] = (
                alldata["LP"]["GenPvar"][datagen].X * baseMVA
            )  # Generator real power injection
            resgen["Qg"] = 0  # Generator reactive power injection is always 0 for DC

        # Branch solution data
        for branchindex in result["branch"]:
            resbranch = result["branch"][branchindex]
            databranch = alldata["branches"][branchindex]
            # Generate new values
            # Real power injected into "from" end of branch and "to" end of branch are the same for DC
            resbranch["Pf"] = resbranch["Pt"] = (
                alldata["LP"]["Pvar_f"][databranch].X * baseMVA
            )
            resbranch["switching"] = 1  # Default value for switching

    # Set AC solution data
    if opftype in [OpfType.AC, OpfType.IV]:
        # Bus solution data only available if we have rectangular formulation or polar
        if alldata["use_ef"]:
            for busindex in result["bus"]:
                resbus = result["bus"][busindex]
                databus = alldata["buses"][busindex]
                # Override old values
                # Voltage magnitude is root of cvar because cvar = square of voltage magnitude given as e^2 + f^2
                if alldata["doiv"]:  # doiv makes sure that e, f variables are present
                    resbus["Vm"] = math.sqrt(
                        alldata["LP"]["evar"][databus].X ** 2
                        + alldata["LP"]["fvar"][databus].X ** 2
                    )
                else:
                    resbus["Vm"] = math.sqrt(alldata["LP"]["cvar"][databus].X)

            if alldata["doiv"]:
                # Need to fill cvar[branch] dictionary to compute angles for IV
                # Note that there is no cvar dictionary for IV!!!
                cvar = {}
                for branch in alldata["branches"].values():
                    fbus = alldata["buses"][branch.f]
                    tbus = alldata["buses"][branch.t]
                    cvar[branch] = (
                        alldata["LP"]["evar"][fbus].X * alldata["LP"]["evar"][tbus].X
                        + alldata["LP"]["fvar"][fbus].X * alldata["LP"]["fvar"][tbus].X
                    )

                alldata["LP"]["cvar"] = cvar

            compute_voltage_angles(alldata, result)

        if alldata["dopolar"] and not alldata["doiv"]:
            for busindex in result["bus"]:
                resbus = result["bus"][busindex]
                databus = alldata["buses"][busindex]
                # Override old values
                resbus["Vm"] = alldata["LP"]["vvar"][databus].X  # Voltage magnitude
                resbus["Va"] = alldata["LP"]["thetavar"][databus].X  # Voltage angle

        # Generator solution data
        for genindex in result["gen"]:
            resgen = result["gen"][genindex]
            datagen = alldata["gens"][genindex]
            # Override old values
            resgen["Pg"] = (
                alldata["LP"]["GenPvar"][datagen].X * baseMVA
            )  # Generator real power injection
            resgen["Qg"] = (
                alldata["LP"]["GenQvar"][datagen].X * baseMVA
            )  # Generator reactive power injection

        # Branch solution data
        for branchindex in result["branch"]:
            resbranch = result["branch"][branchindex]
            databranch = alldata["branches"][branchindex]
            # Generate new values
            # Real power injected into "from" end of branch and "to" end of branch are the same for DC
            resbranch["Pf"] = (
                alldata["LP"]["Pvar_f"][databranch].X * baseMVA
            )  # AC branch real power injected into "from" end of branch
            resbranch["Pt"] = (
                alldata["LP"]["Pvar_t"][databranch].X * baseMVA
            )  # AC branch real power injected into "to" end of branch
            resbranch["Qf"] = (
                alldata["LP"]["Qvar_f"][databranch].X * baseMVA
            )  # AC branch reactive power injected into "from" end of branch
            resbranch["Qt"] = (
                alldata["LP"]["Qvar_t"][databranch].X * baseMVA
            )  # AC branch reactive power injected into "to" end of branch
            resbranch["switching"] = 1  # Default value for switching

    # Set MIP fields
    if alldata["branchswitching_mip"] or (
        opftype == OpfType.AC and alldata["branchswitching_comp"]
    ):
        for branchindex in result["branch"]:
            resbranch = result["branch"][branchindex]
            databranch = alldata["branches"][branchindex]
            resbranch["switching"] = (
                1 if alldata["MIP"]["zvar"][databranch].X > 0.5 else 0
            )


def compute_voltage_angles(alldata, result):
    """
    Helper function to compute voltage angles out of previously computed
    voltage magnitudes for a given AC OPF solution

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param result: Result dictionary which is constructed for the user
    :type result: dict
    """

    # After setting all voltage magnitudes, we can compute the voltage angle
    # Set voltage angle of reference bus to 0 (this is arbitrary)
    buses = alldata["buses"]
    branches = alldata["branches"]
    busindex = alldata["refbus"]
    cvar = alldata["LP"]["cvar"]
    result["bus"][busindex]["Va"] = 0
    # Next, compute the angle of all other buses starting from the reference bus
    # k = "from" m = "to"
    # cvar[km] = (voltage mag at k) * (voltage mag at m) * cos(thetak - thetam)
    busesdone = {busindex}
    busesnext = set()
    for branchindex in buses[busindex].frombranchids.values():
        b = branches[branchindex]
        # Already computed angle is always at first place of tuple
        busesnext.add((b.f, b.t, branchindex, "f"))
    for branchindex in buses[busindex].tobranchids.values():
        b = branches[branchindex]
        # Already computed angle is always at first place of tuple
        busesnext.add((b.t, b.f, branchindex, "t"))

    while len(busesdone) != alldata["numbuses"]:
        # Get a new bus
        nextb = busesnext.pop()
        while nextb[1] in busesdone:
            nextb = busesnext.pop()
        nextbusindex = nextb[1]
        nextbus = buses[nextbusindex]
        knownbusindex = nextb[0]
        buses[knownbusindex]
        if alldata["doiv"]:
            # For IV, we filled the values manually so they are not Gurobi variables
            cvarval = cvar[branches[nextb[2]]]
        else:
            cvarval = cvar[branches[nextb[2]]].X

        tmp = cvarval / (
            result["bus"][nextbusindex]["Vm"] * result["bus"][knownbusindex]["Vm"]
        )
        # Safe guard to avoid slightly being out of [-1,1] for numerical reasons
        tmp = max(-1, min(1, tmp))
        res = math.acos(tmp)
        if nextb[3] == "f":
            res -= result["bus"][knownbusindex]["Va"]
            res *= -1
        else:
            res += result["bus"][knownbusindex]["Va"]
        result["bus"][nextbusindex]["Va"] = res
        busesdone.add(nextbusindex)

        for branchindex in nextbus.frombranchids.values():
            b = branches[branchindex]
            # Only add if bus is not done already
            if b.t not in busesdone:
                busesnext.add((b.f, b.t, branchindex, "f"))
        for branchindex in nextbus.tobranchids.values():
            b = branches[branchindex]
            if b.f not in busesdone:
                busesnext.add((b.t, b.f, branchindex, "t"))


def fill_violations_fields(alldata, opftype, result):
    """
    Extends the result dictionary by additional violation information

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param model: Gurobi model which was recently optimized
    :type model: :class: `gurobipy.Model`
    :param opftype: Type of OPF formulation
    :type opftype: :enum: `OpfType`
    :param result: Dictionary holding all case and violation data
    :type result: dict

    :raises ValueError: Unsupported OPF type
    """
    baseMVA = result["baseMVA"]

    if opftype != OpfType.AC:
        raise ValueError("Only AC model type supported for violations.")

    # Bus violations
    for busindex in result["bus"]:
        resbus = result["bus"][busindex]
        databus = alldata["buses"][busindex]
        resbus["Vmviol"] = (
            alldata["violation"]["Vmagviol"][databus] * baseMVA
        )  # Voltage magnitude violation
        resbus["Pviol"] = (
            alldata["violation"]["IPviol"][databus] * baseMVA
        )  # Real injection violation
        resbus["Qviol"] = (
            alldata["violation"]["IQviol"][databus] * baseMVA
        )  # Reactive injection violation

    # Branch limit violations
    for branchindex in result["branch"]:
        resbranch = result["branch"][branchindex]
        databranch = alldata["branches"][branchindex]
        resbranch["limitviol"] = (
            alldata["violation"]["branchlimit"][databranch] * baseMVA
        )
