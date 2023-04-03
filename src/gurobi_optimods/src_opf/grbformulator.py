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
    lpformulator_iv_strictchecker,
)


def construct_and_solve_model(alldata):
    """Construct OPF model and solve it"""

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

    logging.info("\n%s formulation." % opftype.value)

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

        if alldata["strictcheckvoltagesolution"]:
            # check input solution against formulation
            spitoutvector = True
            # TODO-Dan Why do you need the feascode? It's not used anywhere and the function does not even return anything
            feascode = lpformulator_strictchecker(
                alldata, model, spitoutvector, opftype
            )

        if alldata["doslp_polar"]:
            break_exit("slp_polar formulation")  # TODO-Dan why the break_exit?

        sol_count = lpformulator_optimize(alldata, model, opftype)

        endtime = time.time()
        logging.info(
            "Overall time taken (model construction + optimization): %f s."
            % (endtime - starttime)
        )
        logging.info("Solution count: %d." % (sol_count))

        if sol_count > 0:
            lpformulator_examine_solution(alldata, model, opftype)
            objval = model.ObjVal
            # we are using an ordered dict to maintain the index of variables so we can access them by index when plotting solution
            solution = OrderedDict()
            index = 0
            for v in model.getVars():
                if math.fabs(v.X) > 1e-09:
                    logging.info(v.varname + " = " + str(v.X))
                    solution[v.VarName] = v.X
                else:
                    solution[v.VarName] = 0.0
                index += 1

    return solution, objval

    """
    buses        = alldata['buses']
    branches     = alldata['branches']
    gens         = alldata['gens']
    IDtoCountmap = alldata['IDtoCountmap']
    """


def lpformulator_body(alldata, model, opftype):
    """Call the corresponding model construction method"""
    if opftype == OpfType.AC:
        lpformulator_ac_body(alldata, model)
    elif opftype == OpfType.DC:
        lpformulator_dc_body(alldata, model)
    elif opftype == OpfType.IV:
        lpformulator_iv_body(alldata, model)
    else:
        raise ValueError("Unknown OPF type.")


def lpformulator_strictchecker(alldata, model, spitoutvector, opftype):
    if opftype == OpfType.AC:
        lpformulator_ac_strictchecker(alldata, model, spitoutvector)
    elif opftype == OpfType.DC:
        # lpformulator_dc_strictchecker(alldata, model, spitoutvector)
        pass  # TODO-Dan Is there a reason why there is no strict checker for DC (except that it's linear)
    elif opftype == OpfType.IV:
        lpformulator_iv_strictchecker(alldata, model, spitoutvector)
    else:
        raise ValueError("Unknown OPF type.")


def lpformulator_examine_solution(alldata, model, opftype):
    """Call the corresponding model construction method"""
    if opftype == OpfType.AC:
        lpformulator_ac_examine_solution(alldata, model)
    elif opftype == OpfType.DC:
        lpformulator_dc_examine_solution(alldata, model)
    elif opftype == OpfType.IV:
        lpformulator_iv_examine_solution(alldata, model)
    else:
        raise ValueError("Unknown OPF type.")


def lpformulator_optimize(alldata, model, opftype):
    """Optimize constructed model"""

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

    # TODO-Dan To me it looked like only lpformulator_dc_opt had a different if check for no reason. Can you explain?
    # The check for DC was
    # if alldata["branchswitching_mip"]:
    #  The content of the if-clause is the same among DC,AC
    if (
        alldata["usemipstart"]
        and (alldata["branchswitching_mip"] or alldata["branchswitching_comp"])
        or (alldata["branchswitching_mip"] and opftype == OpfType.DC)
    ):
        # logging level needs to be critical because INFO is disabled
        logging.critical("Using mip start with all branches kept on.")
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
        logging.info("\nModel Status: infeasible or unbounded.\n")
        logging.info("Re-optimizing with DualReductions turned off.\n")
        logging.disable(logging.INFO)
        model.Params.DualReductions = 0
        model.optimize()
        logging.disable(logging.NOTSET)

    if model.status == GRB.INFEASIBLE:
        logging.info("\nModel Status: infeasible.\n")
        logging.info("Computing IIS...")
        logging.disable(logging.INFO)
        model.computeIIS()
        logging.disable(logging.NOTSET)
        logging.info("\nIIS computed, writing IIS to file acopfmodel.ilp.")
        model.write("acopfmodel_iis.ilp")

    elif model.status == GRB.UNBOUNDED:
        logging.info("\nModel Status: unbounded.\n")

    elif model.status == GRB.INTERRUPTED:
        logging.info("\nModel Status: interrupted.\n")

    elif model.status == GRB.OPTIMAL:
        logging.info("\nModel Status: optimal.\n")

    # Only print objective value and solution quality if at least
    # one feasible point is available
    if model.SolCount > 0:
        logging.info("Objective value = %.8e." % model.objVal)
        logging.disable(logging.INFO)
        model.printQuality()
        logging.disable(logging.NOTSET)

    return model.SolCount


def lpformulator_setup(alldata, opftype):
    """Helper function to handle specific settings"""

    # Additional setup is only meant for AC and IV
    if opftype == OpfType.DC:
        return

    logging.info("Auxiliary setup.")

    alldata["maxdispersion_rad"] = (math.pi / 180.0) * alldata["maxdispersion_deg"]

    if alldata["dopolar"]:
        logging.info("  Polar formulation, shutting down incompatible options:")
        alldata["use_ef"] = False
        alldata["useconvexformulation"] = False
        alldata["skipjabr"] = True
        logging.info("    use_ef, useconvexformulation, jabr.")

    if alldata["voltsfilename"] != None:
        grbreadvoltsfile(alldata)

    # The following settings are only for AC
    if opftype != OpfType.AC:
        return

    alldata["maxphasediff_rad"] = (math.pi / 180.0) * alldata["maxphasediff_deg"]

    branches = alldata["branches"]
    numbranches = alldata["numbranches"]
    buses = alldata["buses"]

    if alldata["usemaxphasediff"]:
        logging.info(
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

        logging.info("Updated %d maxangle constraints." % (count))


def writempsfile(alldata, model, filename):
    """Helper function for debugging"""
    logging.info("Writing mpsfile to %s." % (filename))
    model.write(filename)


def writemipstart(alldata):

    filename = "mipstart.mst"
    f = open(filename, "w")
    logging.info("Writing mipstart in file %s." % filename)

    zvar = alldata["LP"]["zvar"]
    branches = alldata["branches"]
    numbranches = alldata["numbranches"]
    for j in range(1, 1 + numbranches):
        branch = branches[j]
        f.write("{} 1.0\n".format(zvar[branch].Varname))

    f.close()
