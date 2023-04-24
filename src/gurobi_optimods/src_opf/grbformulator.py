import math
import time
import logging
import numpy as np
import gurobipy as gp
from gurobipy import GRB

from .utils import OpfType
from .grbfile import grbreadvoltsfile
from .grbformulator_ac import (
    lpformulator_ac_body,
    lpformulator_ac_strictchecker,
)
from .grbformulator_dc import lpformulator_dc_body
from .grbformulator_iv import (
    lpformulator_iv_body,
)


def construct_and_solve_model(alldata):
    """
    Construct OPF model and solve it

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    """

    logger = logging.getLogger("OpfLogger")
    opftype = None
    if alldata["doac"]:
        opftype = OpfType.AC
    elif alldata["dodc"]:
        opftype = OpfType.DC
    elif alldata["doiv"]:
        opftype = OpfType.IV
    else:
        # This should have been checked before already but it can't hurt to have this error
        raise ValueError(
            "Illegal option combination. Have to use exactly 1 of options [doac, dodc, doiv]."
        )

    logger.info("\n%s formulation." % opftype.value)

    starttime = time.time()

    sol_count = 0
    solution = None
    objval = None
    modelname = opftype.value + "_Formulation_Model"

    # Handle special settings
    if opftype != OpfType.DC:
        lpformulator_setup(alldata, opftype)

    # Create model
    with gp.Env() as env, gp.Model(modelname, env=env) as model:
        # Add model variables and constraints
        lpformulator_body(alldata, model, opftype)

        # Update to get correct model stats
        model.update()
        logger.info(
            "Constructed %sOPF model with %d variables and %d constraints.\n"
            % (opftype.value, model.NumVars, model.NumConstrs)
        )

        # Write model to file if requested by user
        if alldata["lpfilename"] != None:
            model.write(alldata["lpfilename"])
            logger.info("Wrote LP to " + alldata["lpfilename"])

        # TODO This will be an own API function
        if alldata["strictcheckvoltagesolution"]:
            # check input solution against formulation
            spitoutvector = True
            # TODO-Dan Why do you need the feascode? It's not used anywhere and the function does not even return anything
            feascode = lpformulator_strictchecker(
                alldata, model, spitoutvector, opftype
            )

        # if alldata["doslp_polar"]: # TODO-Dan Is this work in progress?
        #    break_exit("slp_polar formulation")  # TODO-Dan why the break_exit?

        # Solve the OPF model
        sol_count = lpformulator_optimize(alldata, model, opftype)

        endtime = time.time()
        logger.info(
            "Overall time taken (model construction + optimization): %f s."
            % (endtime - starttime)
        )
        logger.info("Solution count: %d." % (sol_count))

        # Need to turn Gurobi solution into a dictionary following MATPOWER notation
        solution = turn_solution_into_result_dict(alldata, model, opftype)

        if solution["success"] == 1:
            objval = solution["f"]
        else:
            objval = GRB.INFINITY

    return solution, objval

    """
    buses        = alldata['buses']
    branches     = alldata['branches']
    gens         = alldata['gens']
    IDtoCountmap = alldata['IDtoCountmap']
    """


def lpformulator_body(alldata, model, opftype):
    """
    Call the corresponding model construction method

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    model : gurobipy.Model
        Gurobi model to be constructed
    opftype : OpfType
        Type of OPF formulation
    """
    if opftype == OpfType.AC:
        lpformulator_ac_body(alldata, model)
    elif opftype == OpfType.DC:
        lpformulator_dc_body(alldata, model)
    elif opftype == OpfType.IV:
        lpformulator_iv_body(alldata, model)
    else:
        raise ValueError("Unknown OPF type.")


def lpformulator_strictchecker(alldata, model, spitoutvector, opftype):
    """
    Call the corresponding strictchecker function

    TODO-Dan How do we use this function? It says "check input solution against formulation".
             Which input solution? Where does the user provide the input solution? Is it in the voltage file?

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    model : gurobipy.Model
        Gurobi model to be constructed
    spitoutvector: boolean
        # TODO-Dan What does it do?
    opftype : OpfType
        Type of OPF formulation
    """
    if opftype == OpfType.DC:
        # lpformulator_dc_strictchecker(alldata, model, spitoutvector)
        pass  # TODO-Dan Is there a reason why there is no strict checker for DC (except that it's linear)
    elif opftype == OpfType.AC or opftype == OpfType.IV:
        lpformulator_ac_strictchecker(alldata, model, spitoutvector)
    else:
        raise ValueError("Unknown OPF type.")


def lpformulator_optimize(alldata, model, opftype):
    """
    Optimize constructed model

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    model : gurobipy.Model
        Constructed Gurobi model
    opftype : OpfType
        Type of OPF formulation

    Returns
    -------
    Number of found solutions
    """

    logger = logging.getLogger("OpfLogger")
    # Disable logging handler to get Gurobi output
    logging.disable(logging.INFO)
    model.params.LogFile = alldata["logfile"]
    # Specific settings for better convergence

    if opftype != OpfType.DC:
        if alldata["use_ef"]:
            model.params.NonConvex = 2

        # Specific settings for better convergence
        model.Params.MIPGap = 1.0e-3
        model.Params.OptimalityTol = 1.0e-3
    else:
        model.Params.MIPGap = 1.0e-4
        model.Params.OptimalityTol = 1.0e-4

    if alldata["gurobiparamfile"] != None:
        model.read(alldata["gurobiparamfile"])

    # model.Params.SolutionLimit = 20  # TODO-Dan This setting was set for DC. Why do we need a solution limit for DC?

    # TODO-Dan Why is DC handled differently? It looks like the user does not have to set the usemipstart setting but only
    #          the branchswitching_mip setting if they are solving DC
    if (
        alldata["usemipstart"]
        and (alldata["branchswitching_mip"] or alldata["branchswitching_comp"])
        or (alldata["branchswitching_mip"] and opftype == OpfType.DC)
    ):
        # logging level needs to be critical because INFO is disabled
        logger.critical("Using mip start with all branches kept on.")
        # mip start
        zvar = alldata["MIP"]["zvar"]
        branches = alldata["branches"]
        numbranches = alldata["numbranches"]
        for j in range(1, 1 + numbranches):
            branch = branches[j]
            zvar[branch].Start = 1.0
            # Jarek, there is some strange behavior here
            # print(zvar[branch], ' ',zvar[branch].Start)
            # TODO-Dan What strange behavior?

        # writemipstart(alldata)

    model.optimize()

    # Activate logging handlers
    logging.disable(logging.NOTSET)

    # Check model status and re-optimize or try computing an IIS if necessary
    if model.status == GRB.INF_OR_UNBD:
        logger.info("\nModel Status: infeasible or unbounded.\n")
        logger.info("Re-optimizing with DualReductions turned off.\n")
        logging.disable(logging.INFO)
        model.Params.DualReductions = 0
        model.optimize()
        logging.disable(logging.NOTSET)

    if model.status == GRB.INFEASIBLE:
        logger.info("\nModel Status: infeasible.\n")
        logger.info("Computing IIS...")
        logging.disable(logging.INFO)
        model.computeIIS()
        logging.disable(logging.NOTSET)
        logger.info("\nIIS computed, writing IIS to file acopfmodel.ilp.")
        model.write("acopfmodel_iis.ilp")

    if model.status == GRB.NUMERIC:
        logger.info("\nModel Status: Numerically difficult model.\n")
        logger.info(
            "Re-optimizing with settings focused on improving numerical stability.\n"
        )
        logging.disable(logging.INFO)
        model.Params.NumericFocus = 2
        model.Params.BarHomogeneous = 1
        model.reset()
        model.optimize()
        logging.disable(logging.NOTSET)

    elif model.status == GRB.UNBOUNDED:
        logger.info("\nModel Status: unbounded.\n")

    elif model.status == GRB.INTERRUPTED:
        logger.info("\nModel Status: interrupted.\n")

    elif model.status == GRB.OPTIMAL:
        logger.info("\nModel Status: optimal.\n")

    # Only print objective value and solution quality if at least
    # one feasible point is available
    if model.SolCount > 0:
        logger.info("Objective value = %.8e." % model.objVal)
        logging.disable(logging.INFO)
        model.printQuality()
        logging.disable(logging.NOTSET)

    return model.SolCount


def lpformulator_setup(alldata, opftype):
    """
    Helper function to handle specific settings before starting optimization

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    opftype : OpfType
        Type of OPF formulation
    """

    logger = logging.getLogger("OpfLogger")
    # Additional setup is only meant for AC and IV
    if opftype == OpfType.DC:
        return

    logger.info("Auxiliary setup.")

    alldata["maxdispersion_rad"] = (math.pi / 180.0) * alldata["maxdispersion_deg"]

    if alldata["dopolar"]:
        logger.info("  Polar formulation, shutting down incompatible options:")
        alldata["use_ef"] = False
        alldata["useconvexformulation"] = False
        alldata["skipjabr"] = True
        logger.info("    use_ef, useconvexformulation, jabr.")

    if alldata["voltsfilename"] != None:
        grbreadvoltsfile(alldata)

    if alldata["doiv"] and not alldata["use_ef"]:
        alldata["use_ef"] = True
        logger.info("  IV formulation requires use_ef. Turning it on.")

    # The following settings are only for AC
    if opftype != OpfType.AC:
        return

    alldata["maxphasediff_rad"] = (math.pi / 180.0) * alldata["maxphasediff_deg"]

    branches = alldata["branches"]
    numbranches = alldata["numbranches"]
    buses = alldata["buses"]

    if alldata["usemaxphasediff"]:
        logger.info(
            "Applying max phase diff of %g degrees." % (alldata["maxphasediff_deg"])
        )

        maxrad = alldata["maxphasediff_rad"]

        count = 0

        for j in range(1, 1 + numbranches):
            branch = branches[j]

            if branch.maxangle_rad > maxrad:
                branch.maxangle_rad = maxrad
                count += 1
            if branch.minangle_rad < -maxrad:
                branch.minangle_rad = -maxrad

        logger.info("Updated %d maxangle constraints." % (count))


def turn_solution_into_result_dict(alldata, model, opftype):
    """
    Function to turn a Gurobi solution into an OPF dictionary in MATPOWER notation.
    If no solution is present, the result dictionary "success" value is 0.

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    model : gurobipy.Model
        Gurobi model which was recently optimized
    """

    # Reconstruct case data from our data
    result = {}
    result["baseMVA"] = alldata["baseMVA"]
    baseMVA = result["baseMVA"]
    result["bus"] = {}
    result["gen"] = {}
    result["branch"] = {}
    result["gencost"] = {}

    # buses
    buses = alldata["buses"]
    index = 1
    for b in buses.values():
        # bus_i type    Pd  Qd  Gs  Bs  area    Vm  Va  baseKV  zone    Vmax    Vmin
        # we don't use area, Vm, Va, zone but we save them for consistency with MATPOWER
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
            # "mu": # balance constraints
        }
        result["bus"][index] = matbus
        index += 1

    # generators and gen costs
    gens = alldata["gens"]
    index = 1
    for g in gens.values():
        # generator data
        # bus   Pg  Qg  Qmax    Qmin    Vg  mBase   status  Pmax    Pmin    Pc1 Pc2 Qc1min  Qc1max  Qc2min  Qc2max  ramp_agc    ramp_10 ramp_30 ramp_q  apf
        # we don't use Vg, mBase, and everything after Pmin but we save them for consistency with MATPOWER
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

        # generator cost data
        #   1   startup shutdown    n   x1  y1  ... xn  yn
        #   2   startup shutdown    n   c(n-1)  ... c0
        # we don't use startup and shutdown but we save them for consistency with MATPOWER
        matgencost = {
            "costtype": g.costtype,
            "startup": g.startup,
            "shutdown": g.shutdown,
            "n": g.costdegree + 1,
            "costvector": gencost,
        }
        result["gencost"][index] = matgencost
        index += 1

    # branches
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
            # "switching": 1 alldata["MIP"]["zholder"] index is number of branch
        }
        result["branch"][index] = matbranch
        index += 1

    # Now extend the result dictionary by additional information derived from Gurobi solution
    # See MATPOWER-manual for more details https://matpower.org/docs/MATPOWER-manual.pdf
    result["et"] = model.Runtime
    if model.SolCount == 0:
        # We did not find any feasible solution
        result["success"] = 0
        return result

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
                ].Pi  # shadow prices of balance constraints

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
        if (
            alldata["use_ef"]
            and not alldata["useconvexformulation"]
            and not alldata["doiv"]
        ):  # TODO check how to compute voltage magnitude and angle values for IV
            for busindex in result["bus"]:
                resbus = result["bus"][busindex]
                databus = alldata["buses"][busindex]
                # Override old values
                # Voltage magnitude is root of cvar because cvar = square of voltage magnitude given as e^2 + f^2
                resbus["Vm"] = math.sqrt(alldata["LP"]["cvar"][databus].X)

            compute_voltage_angles(alldata, result)

        if alldata["dopolar"]:
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
    if alldata["branchswitching_mip"]:
        for branchindex in result["branch"]:
            resbranch = result["branch"][branchindex]
            databranch = alldata["branches"][branchindex]
            resbranch["switching"] = (
                1 if alldata["MIP"]["zvar"][databranch].X > 0.5 else 0
            )

    return result


def compute_voltage_angles(alldata, result):
    """
    Helper function to compute voltage angles out of previously computed
    voltage magnitudes for a given AC OPF solution

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    result : dictionary
        Result dictionary which is constructed for the user
    """

    # After setting all voltage magnitudes, we can compute the voltage angle
    # Set voltage angle of reference bus to 0 (this is arbitrary)
    buses = alldata["buses"]
    branches = alldata["branches"]
    busindex = alldata["refbus"]
    result["bus"][busindex]["Va"] = 0
    # Next, compute the angle of all other buses starting from the reference bus
    # k = "from" m = "to"
    # cvar[km] = (voltage mag at k) * (voltage mag at m) * cos(thetak - thetam)
    busesdone = {busindex}
    busesnext = set()
    for branchindex in buses[busindex].frombranchids.values():
        b = branches[branchindex]
        # already computed angle is always at first place of tuple
        busesnext.add((b.f, b.t, branchindex, "f"))
    for branchindex in buses[busindex].tobranchids.values():
        b = branches[branchindex]
        # already computed angle is always at first place of tuple
        busesnext.add((b.t, b.f, branchindex, "t"))

    while len(busesdone) != alldata["numbuses"]:
        # Get a new bus
        next = busesnext.pop()
        while next[1] in busesdone:
            next = busesnext.pop()
        nextbusindex = next[1]
        nextbus = buses[nextbusindex]
        knownbusindex = next[0]
        knownbus = buses[knownbusindex]
        cvarval = alldata["LP"]["cvar"][branches[next[2]]].X
        res = math.acos(
            cvarval
            / (result["bus"][nextbusindex]["Vm"] * result["bus"][knownbusindex]["Vm"])
        )
        if next[3] == "f":
            res -= result["bus"][knownbusindex]["Va"]
            res *= -1
        else:
            res += result["bus"][knownbusindex]["Va"]
        result["bus"][nextbusindex]["Va"] = res
        # print("setting voltage angle of bus %d to value %f"%(nextbusindex,res))
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


def writemipstart(alldata):
    """
    Helper function for writing a mip start to mipstart.mst file.
    Mainly used for debugging numerically difficult cases

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data

    Returns
    -------
    Number of found solutions
    """

    logger = logging.getLogger("OpfLogger")
    filename = "mipstart.mst"
    f = open(filename, "w")
    logger.info("Writing mipstart in file %s." % filename)

    zvar = alldata["MIP"]["zvar"]
    branches = alldata["branches"]
    numbranches = alldata["numbranches"]
    for j in range(1, 1 + numbranches):
        branch = branches[j]
        f.write("{} 1.0\n".format(zvar[branch].Varname))

    f.close()
