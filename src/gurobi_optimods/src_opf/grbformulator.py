import math
import time
import logging
import numpy as np
import gurobipy as gp
from collections import OrderedDict
from gurobipy import GRB

from .utils import OpfType
from .grbfile import grbreadvoltsfile
from .grbformulator_ac import (
    lpformulator_ac_body,
    lpformulator_ac_examine_solution,
    lpformulator_ac_strictchecker,
)
from .grbformulator_dc import lpformulator_dc_body, lpformulator_dc_examine_solution
from .grbformulator_iv import (
    lpformulator_iv_body,
    lpformulator_iv_examine_solution,
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

        # TODO This will be an own API function
        if alldata["strictcheckvoltagesolution"]:
            # check input solution against formulation
            spitoutvector = True
            # TODO-Dan Why do you need the feascode? It's not used anywhere and the function does not even return anything
            feascode = lpformulator_strictchecker(
                alldata, model, spitoutvector, opftype
            )

        # if alldata["doslp_polar"]:
        #    break_exit("slp_polar formulation")  # TODO-Dan why the break_exit?

        sol_count = lpformulator_optimize(alldata, model, opftype)

        endtime = time.time()
        logger.info(
            "Overall time taken (model construction + optimization): %f s."
            % (endtime - starttime)
        )
        logger.info("Solution count: %d." % (sol_count))

        if sol_count > 0:
            lpformulator_examine_solution(alldata, model, opftype)
            objval = model.ObjVal
            # we are using an ordered dict to maintain the index of variables so we can access them by index when plotting solution
            solution = OrderedDict()
            index = 0
            for v in model.getVars():
                if math.fabs(v.X) > 1e-09:
                    if (
                        model.NumVars < 60
                    ):  # Only print the solution if it's not too big
                        logger.info(v.varname + " = " + str(v.X))
                    solution[v.VarName] = v.X
                else:
                    solution[v.VarName] = 0.0
                index += 1

            turn_solution_into_mpc_dict(alldata, model)

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


def lpformulator_examine_solution(alldata, model, opftype):
    """
    Call the corresponding solution examination function

    # TODO-Dan I added a comment in grbformulator_ac.py:lpformulator_ac_examine_solution

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    model : gurobipy.Model
        Constructed Gurobi model
    opftype : OpfType
        Type of OPF formulation
    """
    if opftype == OpfType.AC:
        lpformulator_ac_examine_solution(alldata, model)
    elif opftype == OpfType.DC:
        lpformulator_dc_examine_solution(alldata, model)
    elif opftype == OpfType.IV:
        lpformulator_iv_examine_solution(alldata, model)
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

    model.Params.FeasibilityTol = (
        1.0e-6  # TODO-Dan it's Gurobi's default setting, we can remove this
    )

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
        zvar = alldata["LP"]["zvar"]
        branches = alldata["branches"]
        numbranches = alldata["numbranches"]
        for j in range(1, 1 + numbranches):
            branch = branches[j]
            zvar[branch].Start = 1.0
            # Jarek, there is some strange behavior here
            # print(zvar[branch], ' ',zvar[branch].Start)
            # TODO-Dan What strange behavior?

        # writemipstart(alldata)

        zholder = np.zeros(numbranches)
        alldata["MIP"]["zholder"] = zholder
        alldata["MIP"]["solutionfound"] = False
        alldata["MIP"]["bestsolval"] = 1e50
        alldata["MIP"]["solcount"] = 0
        gholder = np.zeros(alldata["numgens"])
        alldata["MIP"]["gholder"] = gholder

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


def turn_solution_into_mpc_dict(alldata, model):

    mpc = {}
    mpc["baseMVA"] = alldata["baseMVA"]
    baseMVA = mpc["baseMVA"]
    mpc["bus"] = {}
    mpc["gen"] = {}
    mpc["branch"] = {}
    mpc["gencost"] = {}
    vars = model.getVars()

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
        }
        print("bus %d P var value: %f" % (b.nodeID, b.Pinjvarind))
        mpc["bus"][index] = matbus
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
        mpc["gen"][index] = matgen
        print("gen %d P var value: %f" % (g.nodeID, g.Pvarind))
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
        mpc["gencost"][index] = matgencost
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
        }
        mpc["branch"][index] = matbranch
        index += 1

    # print(mpc)


def writempsfile(alldata, model, filename):
    """
    Helper function for writing Gurobi models.
    Mainly used for debugging

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    model : gurobipy.Model
        Constructed Gurobi model
    filename : string
        Name of model file

    Returns
    -------
    Number of found solutions
    """
    logger = logging.getLogger("OpfLogger")
    logger.info("Writing model to %s." % (filename))
    model.write(filename)


def writemipstart(alldata):
    """
    Helper function for writing a mip start to mipstart.mst file.
    Mainly used for debugging

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

    zvar = alldata["LP"]["zvar"]
    branches = alldata["branches"]
    numbranches = alldata["numbranches"]
    for j in range(1, 1 + numbranches):
        branch = branches[j]
        f.write("{} 1.0\n".format(zvar[branch].Varname))

    f.close()
