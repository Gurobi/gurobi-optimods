import math
import time
import logging
import numpy as np
import gurobipy as gp
from gurobipy import GRB

from .myutils import break_exit
from .grbfile import grbreadvoltsfile
from .grbgraphical import grbgraphical


def lpformulator_ac(alldata):
    """Formulate ACOPF model and solve it"""

    # TODO-Dan The lpformulator_XX functions all look the same. Couldn't we just use one and have a switch-case over alldata["dodc"]/alldata["doac"]/alldata["doiv"]
    #         Or is there some difference?

    logging.info("\nAC formulation.")

    starttime = time.time()

    # Handle special settings
    lpformulator_setup(alldata)
    sol_count = 0
    solution = None
    objval = None

    # Create model
    with gp.Env() as env, gp.Model("ac_formulation_model", env=env) as model:
        # Add model variables and constraints
        lpformulator_ac_body(alldata, model)

        if alldata["strictcheckvoltagesolution"]:
            # check input solution against formulation
            spitoutvector = True
            # TODO-Dan Why do you need the feascode? It's not used anywhere and the function does not even return anything
            feascode = lpformulator_ac_strictchecker(alldata, model, spitoutvector)

        if alldata["doslp_polar"]:
            break_exit("slp_polar formulation")  # TODO-Dan why the break_exit?

        sol_count = lpformulator_ac_opt(alldata, model)

        endtime = time.time()
        logging.info(
            "Overall time taken (model construction + optimization): %f s."
            % (endtime - starttime)
        )
        logging.info("Solution count: %d." % (sol_count))

        if sol_count > 0:
            lpformulator_ac_examine_solution(alldata, model)
            objval = model.ObjVal
            solution = {}
            for v in model.getVars():
                if math.fabs(v.x) > 1e-09:
                    logging.info(v.varname + " = " + str(v.x))
                    solution[v.VarName] = v.X
                else:
                    solution[v.VarName] = 0

    return solution, objval

    """
    buses        = alldata['buses']
    branches     = alldata['branches']
    gens         = alldata['gens']
    IDtoCountmap = alldata['IDtoCountmap']
    """


def lpformulator_ac_opt(alldata, model):
    """Optimize constructed model"""

    # Disable logging handler to get Gurobi output
    logging.disable(logging.INFO)
    model.params.LogFile = alldata["logfile"]
    # Specific settings for better convergence
    if alldata["use_ef"]:
        model.params.NonConvex = 2

    model.Params.MIPGap = 1.0e-3
    model.Params.OptimalityTol = 1.0e-3
    model.Params.FeasibilityTol = 1.0e-6

    feastol = model.Params.FeasibilityTol
    opttol = model.Params.OptimalityTol
    mipgap = model.Params.MIPGap

    if alldata["usemipstart"] and (
        alldata["branchswitching_mip"] or alldata["branchswitching_comp"]
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

        writemipstart(alldata)

        zholder = np.zeros(numbranches)
        alldata["MIP"]["zholder"] = zholder
        alldata["MIP"]["solutionfound"] = False
        alldata["MIP"]["bestsolval"] = 1e50
        alldata["MIP"]["solcount"] = 0
        gholder = np.zeros(alldata["numgens"])
        alldata["MIP"]["gholder"] = gholder

        # Optimize
        model._vars = model.getVars()
        model.optimize()
    else:
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


def lpformulator_ac_examine_solution(alldata, model):
    """
    TODO-Dan What is this function and how should it be used?
    Currently it does nothing
    """
    # first version -- crude -- but first, some debugging

    branches = alldata["branches"]
    numbranches = alldata["numbranches"]
    IDtoCountmap = alldata["IDtoCountmap"]
    cvar = alldata["LP"]["cvar"]
    buses = alldata["buses"]
    Pvar_f = alldata["LP"]["Pvar_f"]
    Pvar_t = alldata["LP"]["Pvar_t"]
    Qvar_f = alldata["LP"]["Qvar_f"]
    Qvar_t = alldata["LP"]["Qvar_t"]
    minaloss = 100.0
    minalossind = -1
    """
        for j in range(1,1+numbranches):
            branch     = branches[j]
            if branch.isacline:
                f          = branch.f
                t          = branch.t
                count_of_f = IDtoCountmap[f]
                count_of_t = IDtoCountmap[t]
                busf       = buses[count_of_f]
                bust       = buses[count_of_t]
                aLoss = Pvar_f[branch].x + Pvar_t[branch].x
                rLoss = Qvar_f[branch].x + Qvar_t[branch].x

                Sn2 = (Pvar_f[branch].x)*(Pvar_f[branch].x)  +  (Qvar_f[branch].x)*(Qvar_f[branch].x)
                impliedI2 = Sn2/cvar[busf].x
                impliedaLoss = (branch.r)*impliedI2
                actualI2 = branch.ynorm2*(cvar[busf].x + cvar[bust].x - 2*cvar[branch].x)
                actualLoss_r = (branch.r)*actualI2


                logging.info('j %d f %d Pft %8.4f Qft %8.4f cff %8.4f aL %.2e r %.2e iaL %.3e I2 %.3e aL_r %.3e.'%(j,f,Pvar_f[branch].x, Qvar_f[branch].x, cvar[busf].x, aLoss, branch.r, impliedaLoss, actualI2, actualLoss_r))
                #logging.info('          Qft %8.4f Qtf %8.4f rL %.2e x %.2e.'%(Qvar_f[branch].x, Qvar_t[branch].x, rLoss, branch.x))
                if aLoss < minaloss:
                    minaloss = aLoss
                    minalossind = j



        logging.info('minaloss %.3e at %d.'%(minaloss,minalossind))

        """

    logging.info("Done examining solution.\n")


def lpformulator_setup(alldata):
    """Helper function to handle specific settings"""

    logging.info("Auxiliary setup.")

    alldata["maxdispersion_rad"] = (math.pi / 180.0) * alldata["maxdispersion_deg"]
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

    if alldata["dopolar"]:
        logging.info("  Polar formulation, shutting down incompatible options:")
        alldata["use_ef"] = False
        alldata["useconvexformulation"] = False
        alldata["skipjabr"] = True
        logging.info("    use_ef, useconvexformulation, jabr.")

    if alldata["voltsfilename"] != None:
        grbreadvoltsfile(alldata)


def computebalbounds(alldata, bus):
    """Computes active and reactive max and min bus flow balance values"""

    # First let's get max/min generations
    # TODO-Dan What about all those "loud" settings? can we remove them?
    loud = 0  # See below

    baseMVA = alldata["baseMVA"]
    gens = alldata["gens"]
    Pubound = Plbound = 0
    Qubound = Qlbound = 0

    for gencounter in bus.genidsbycount:
        if gens[gencounter].status:
            Pubound += gens[gencounter].Pmax
            Plbound += gens[gencounter].Pmin

            upper = gens[gencounter].Qmax
            lower = gens[gencounter].Qmin
            if bus.nodetype == 3:
                lower = -GRB.INFINITY
                upper = GRB.INFINITY

            Qubound += upper
            Qlbound += lower

        if (
            loud
        ):  # it is worth keeping because a user may want to debug the generator limits that they imposed, which could make a problem infeasible
            # logging.info(" Pubound for %d %f genc %d."%(bus.nodeID, Pubound, gencounter))
            # logging.info(" Plbound for %d %f genc %d."%(bus.nodeID, Plbound, gencounter))
            # TODO-Dan How is this worth keeping if loud cannot be set by the user? The user would have to hack the code to get the debug output
            logging.info(
                " Qubound for %d %f genc %d." % (bus.nodeID, Qubound, gencounter)
            )
            logging.info(
                " Qlbound for %d %f genc %d." % (bus.nodeID, Qlbound, gencounter)
            )

    Pubound -= bus.Pd
    Plbound -= bus.Pd
    Qubound -= bus.Qd
    Qlbound -= bus.Qd

    if bus.nodetype == 4:
        Pubound = Plbound = Qubound = Qlbound = 0

    if loud:  # same as above
        # TODO-Dan same as above
        # logging.info(" Pubound for %d final %f."%(bus.nodeID, Pubound))
        # logging.info(" (Pd was %g)"%bus.Pd)
        # logging.info(" Plbound for %d final %f."%(bus.nodeID, Plbound))
        logging.info(" Qubound for %d final %f." % (bus.nodeID, Qubound))
        logging.info(" (Qd was %g)" % bus.Qd)
        logging.info(" Qlbound for %d final %f." % (bus.nodeID, Qlbound))

    return Pubound, Plbound, Qubound, Qlbound


def lpformulator_ac_body(alldata, model):
    """Helper function for adding variables and constraints to the model"""

    # Create model variables
    lpformulator_ac_create_vars(alldata, model)
    # Create model constraints
    lpformulator_ac_create_constraints(alldata, model)

    model.update()  # Update to get correct model stats
    logging.info(
        "Constructed ACOPF model with %d variables and %d constraints.\n"
        % (model.NumVars, model.NumConstrs)
    )

    if alldata["lpfilename"] != None:
        model.write(
            alldata["lpfilename"]
        )  # FIXME remove.  Jarek: I am using this for debugging, for now
        # TODO-Dan don't forget to remove at some point!
        logging.info("Wrote LP to " + alldata["lpfilename"])

    alldata["model"] = model


def lpformulator_ac_create_vars(alldata, model):
    """Create model variables for ACOPF"""

    logging.info("Creating variables.")

    fixtolerance = 1e-05

    if alldata["fixtolerance"] > 0:
        fixtolerance = alldata["fixtolerance"]

    numbuses = alldata["numbuses"]
    buses = alldata["buses"]
    IDtoCountmap = alldata["IDtoCountmap"]

    # Bus related variables
    cvar = {}
    svar = {}
    Pinjvar = {}
    Qinjvar = {}
    GenPvar = {}
    GenQvar = {}
    gens = alldata["gens"]
    varcount = 0

    for j in range(1, numbuses + 1):
        bus = buses[j]
        # First, injection variables
        maxprod = bus.Vmax * bus.Vmax
        minprod = bus.Vmin * bus.Vmin
        ubound = maxprod
        lbound = minprod

        if alldata["fixcs"] and bus.inputvoltage:
            lbound = bus.inputV * bus.inputV - fixtolerance
            ubound = bus.inputV * bus.inputV + fixtolerance

        cvar[bus] = model.addVar(
            obj=0.0, lb=lbound, ub=ubound, name="c_%d_%d" % (bus.nodeID, bus.nodeID)
        )

        bus.cffvarind = varcount
        varcount += 1

        # csdefslacks to be done
        Plbound = Qlbound = -GRB.INFINITY
        Pubound = Qubound = GRB.INFINITY
        Pubound, Plbound, Qubound, Qlbound = computebalbounds(alldata, bus)

        Pinjvar[bus] = model.addVar(
            obj=0.0, lb=Plbound, ub=Pubound, name="IP_%d" % bus.nodeID
        )
        bus.Pinjvarind = varcount
        # comment: Pinjvar is the variable modeling total active power injected by bus j into the branches incident with j
        varcount += 1

        Qinjvar[bus] = model.addVar(
            obj=0.0, lb=Qlbound, ub=Qubound, name="IQ_%d" % bus.nodeID
        )
        bus.Qinjvarind = varcount
        # TODO-Dan Why is the varcount += 1 in the below comment? Should it be uncommented?
        # comment: Qinjvar is the variable modeling total reactive power injected by bus j into the branches incident with j           varcount += 1

        # Next, generator variables
        for genid in bus.genidsbycount:
            gen = gens[genid]
            lower = gen.Pmin * gen.status
            upper = gen.Pmax * gen.status
            # if bus.nodetype == 3:
            #  upper = GRB.INFINITY
            #  lower = -GRB.INFINITY  #ignoring slack bus
            GenPvar[gen] = model.addVar(
                obj=0.0, lb=lower, ub=upper, name="GP_%d_%d" % (gen.count, gen.nodeID)
            )
            gen.Pvarind = varcount
            varcount += 1

            lower = gen.Qmin * gen.status
            upper = gen.Qmax * gen.status
            if bus.nodetype == 3:
                lower = -GRB.INFINITY
                upper = GRB.INFINITY

            GenQvar[gen] = model.addVar(
                obj=0.0, lb=lower, ub=upper, name="GQ_%d_%d" % (gen.count, gen.nodeID)
            )
            gen.Qvarind = varcount
            varcount += 1

    # Branch related variables
    branches = alldata["branches"]
    numbranches = alldata["numbranches"]
    for j in range(1, 1 + numbranches):
        branch = branches[j]
        f = branch.f
        t = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]
        busf = buses[count_of_f]
        bust = buses[count_of_t]
        maxprod = buses[count_of_f].Vmax * buses[count_of_t].Vmax
        minprod = buses[count_of_f].Vmin * buses[count_of_t].Vmin

        # Assumption 1.  zero angle difference is always allowed! More precisely minangle_rad <= 0 and maxaxangle_rad >= 0
        if branch.maxangle_rad < 0 or branch.minangle_rad > 0:
            logging.error(
                " --- Broken assumption 1: branch j %d f %d t %d minanglerad %f maxanglerad %f."
                % (j, f, t, branch.minangle_rad, branch.maxangle_rad)
            )
            raise ValueError(
                "Broken assumption 1: branch j %d f %d t %d minanglerad %f maxanglerad %f."
                % (j, f, t, branch.minangle_rad, branch.maxangle_rad)
            )

        ubound = ubasic = maxprod
        lbound = lbasic = -maxprod
        maxanglerad = branch.maxangle_rad
        minanglerad = branch.minangle_rad

        # Cosine
        if maxanglerad <= 0.5 * math.pi:
            # In this case minangle <= 0
            if minanglerad >= -0.5 * math.pi:
                lbound = minprod * min(math.cos(maxanglerad), math.cos(minanglerad))
            elif minanglerad >= -math.pi:
                lbound = maxprod * math.cos(minangle_rad)  # Which is negative
            elif minanglerad >= -1.5 * math.pi:
                lbound = -maxprod
            else:
                lbound = -maxprod

        elif maxanglerad <= math.pi:
            if minanglerad >= -0.5 * math.pi:
                lbound = maxprod * math.cos(maxanglerad)
            elif minanglerad >= -math.pi:
                lbound = maxprod * min(math.cos(maxanglerad), math.cos(minanglerad))
            elif minanglerad >= -1.5 * math.pi:
                lbound = -maxprod
            else:
                lbound = -maxprod

        elif maxanglerad <= 1.5 * math.pi:
            lbound = -maxprod

        elif maxanglerad <= 2 * math.pi:
            lbound = -maxprod

        else:
            ubound = maxprod
            lbound = -maxprod

        if branch.inputcs:
            ubound = branch.inputc + fixtolerance
            lbound = branch.inputc - fixtolerance

        cvar[branch] = model.addVar(
            obj=0.0,
            lb=lbound,
            ub=ubound,
            name="c_%d_%d_%d" % (j, busf.nodeID, bust.nodeID),
        )
        branch.cftvarind = varcount
        varcount += 1

        # Sine
        if maxanglerad <= math.pi / 2:
            ubound = maxprod * math.sin(maxanglerad)

            if minanglerad >= -0.5 * math.pi:
                lbound = maxprod * math.sin(minanglerad)
            elif minanglerad >= -math.pi:
                lbound = -maxprod
            elif minanglerad >= -1.5 * math.pi:
                ubound = maxprod * max(math.sin(maxanglerad), math.sin(minanglerad))
                lbound = -maxprod
            else:
                ubound = maxprod
                lbound = -maxprod

        elif maxanglerad <= math.pi:
            ubound = maxprod

            if minanglerad >= -0.5 * math.pi:
                lbound = maxprod * math.sin(minanglerad)
            elif minanglerad >= -math.pi:
                lbound = -maxprod
            elif minanglerad >= -1.5 * math.pi:
                lbound = -maxprod
            else:
                lbound = -maxprod

        elif maxanglerad <= 1.5 * math.pi:
            ubound = maxprod

            if minanglerad >= -0.5 * math.pi:
                lbound = maxprod * min(math.sin(maxanglerad), math.sin(minanglerad))
            elif minanglerad >= -math.pi:
                lbound = -maxprod
            elif minanglerad >= -1.5 * math.pi:
                lbound = -maxprod
            else:
                lbound = -maxprod
        else:
            ubound = maxprod
            lbound = -maxprod

        svar[branch] = model.addVar(
            obj=0.0,
            lb=lbound,
            ub=ubound,
            name="s_%d_%d_%d" % (j, busf.nodeID, bust.nodeID),
        )
        branch.sftvarind = varcount
        varcount += 1

    Pvar_f = {}
    Qvar_f = {}
    Pvar_t = {}
    Qvar_t = {}

    twinPvar_f = {}
    twinPvar_t = {}
    twinQvar_f = {}
    twinQvar_t = {}

    for j in range(1, 1 + numbranches):
        branch = branches[j]
        f = branch.f
        t = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]
        busf = buses[count_of_f]
        bust = buses[count_of_t]
        if branch.constrainedflow:
            ubound = branch.limit
        else:
            ubound = 2 * (
                abs(alldata["summaxgenP"]) + abs(alldata["summaxgenQ"])
            )  # Generous: assumes line charging up to 100%. However it still amounts to an assumption.
        lbound = -ubound

        Pvar_f[branch] = model.addVar(
            obj=0.0,
            lb=lbound,
            ub=ubound,
            name="P_%d_%d_%d" % (j, busf.nodeID, bust.nodeID),
        )
        branch.Pftvarind = varcount
        varcount += 1

        if alldata["branchswitching_mip"]:
            twinPvar_f[branch] = model.addVar(
                obj=0.0,
                lb=lbound,
                ub=ubound,
                name="twinP_%d_%d_%d" % (j, busf.nodeID, bust.nodeID),
            )
            branch.twinPftvarind = varcount
            varcount += 1

        Pvar_t[branch] = model.addVar(
            obj=0.0,
            lb=lbound,
            ub=ubound,
            name="P_%d_%d_%d" % (j, bust.nodeID, busf.nodeID),
        )
        branch.Ptfvarind = varcount
        varcount += 1

        if alldata["branchswitching_mip"]:
            twinPvar_t[branch] = model.addVar(
                obj=0.0,
                lb=lbound,
                ub=ubound,
                name="twinP_%d_%d_%d" % (j, bust.nodeID, busf.nodeID),
            )
            branch.twinPtfvarind = varcount
            varcount += 1

        Qvar_f[branch] = model.addVar(
            obj=0.0,
            lb=lbound,
            ub=ubound,
            name="Q_%d_%d_%d" % (j, busf.nodeID, bust.nodeID),
        )
        branch.Qftvarind = varcount
        varcount += 1

        if alldata["branchswitching_mip"]:
            twinQvar_f[branch] = model.addVar(
                obj=0.0,
                lb=lbound,
                ub=ubound,
                name="twinQ_%d_%d_%d" % (j, busf.nodeID, bust.nodeID),
            )
            branch.twinQftvarind = varcount
            varcount += 1

        Qvar_t[branch] = model.addVar(
            obj=0.0,
            lb=lbound,
            ub=ubound,
            name="Q_%d_%d_%d" % (j, bust.nodeID, busf.nodeID),
        )
        branch.Qtfvarind = varcount
        varcount += 1

        if alldata["branchswitching_mip"]:
            twinQvar_t[branch] = model.addVar(
                obj=0.0,
                lb=lbound,
                ub=ubound,
                name="twinQ_%d_%d_%d" % (j, bust.nodeID, busf.nodeID),
            )
            branch.twinQtfvarind = varcount
            varcount += 1

    zvar = {}
    if alldata["branchswitching_mip"] or alldata["branchswitching_comp"]:
        logging.info("Adding branch switching variables.")
        for j in range(1, 1 + numbranches):
            branch = branches[j]
            f = branch.f
            t = branch.t
            lbound = 0
            ubound = 1
            zvar[branch] = model.addVar(
                lb=lbound,
                ub=ubound,
                obj=0.0,
                vtype=GRB.INTEGER,
                name="z_%d_%d_%d" % (j, f, t),
            )
            branch.switchvarind = varcount
            varcount += 1
    alldata["LP"]["zvar"] = zvar

    lincostvar = model.addVar(
        obj=1.0, lb=-GRB.INFINITY, ub=GRB.INFINITY, name="lincost"
    )
    alldata["LP"]["lincostvar"] = lincostvar
    alldata["LP"]["lincostvarind"] = varcount
    varcount += 1

    if alldata["usequadcostvar"]:
        quadcostvar = model.addVar(obj=1.0, lb=0, ub=GRB.INFINITY, name="quadcost")
        alldata["LP"]["quadcostvar"] = quadcostvar
        alldata["LP"]["quadcostvarind"] = varcount
        varcount += 1

    constobjval = 0
    for gen in gens.values():
        if gen.status > 0:
            constobjval += gen.costvector[gen.costdegree]

    constvar = model.addVar(obj=constobjval, lb=1.0, ub=1.0, name="constant")
    alldata["LP"]["constvar"] = constvar
    varcount += 1

    # Powerflow variables
    if alldata["use_ef"] and alldata["useconvexformulation"] == False:
        lpformulator_ac_create_efvars(alldata, model, varcount)

    if alldata["dopolar"]:
        lpformulator_ac_create_polar_vars(alldata, model, varcount)

    # Save variable data
    alldata["LP"]["cvar"] = cvar
    alldata["LP"]["svar"] = svar
    alldata["LP"]["Pvar_f"] = Pvar_f
    alldata["LP"]["Pvar_t"] = Pvar_t
    alldata["LP"]["Qvar_f"] = Qvar_f
    alldata["LP"]["Qvar_t"] = Qvar_t
    if alldata["branchswitching_mip"]:
        alldata["LP"]["twinPvar_f"] = twinPvar_f
        alldata["LP"]["twinPvar_t"] = twinPvar_t
        alldata["LP"]["twinQvar_f"] = twinQvar_f
        alldata["LP"]["twinQvar_t"] = twinQvar_t

    alldata["LP"]["GenPvar"] = GenPvar
    alldata["LP"]["GenQvar"] = GenQvar
    alldata["LP"]["Pinjvar"] = Pinjvar
    alldata["LP"]["Qinjvar"] = Qinjvar


def lpformulator_ac_create_polar_vars(alldata, model, varcount):
    """Create polar variables for ACOPF"""

    vvar = {}
    thetavar = {}
    numbuses = alldata["numbuses"]
    buses = alldata["buses"]
    IDtoCountmap = alldata["IDtoCountmap"]
    numbranches = alldata["numbranches"]
    branches = alldata["branches"]
    newvarcount = 0
    fixtolerance = alldata["fixtolerance"]

    logging.info("  Creating variables for polar formulation.")

    for j in range(1, numbuses + 1):
        bus = buses[j]
        ubound = bus.Vmax
        lbound = bus.Vmin

        if bus.inputvoltage:
            candidatelbound = bus.inputV - fixtolerance
            candidateubound = bus.inputV + fixtolerance

            lbound = max(lbound, candidatelbound)
            ubound = min(ubound, candidateubound)

        vvar[bus] = model.addVar(
            obj=0.0, lb=lbound, ub=ubound, name="v_" + str(bus.nodeID)
        )
        bus.vvarind = varcount + newvarcount

        ubound = 2 * math.pi
        lbound = -ubound

        if bus.inputvoltage:
            candidatelbound = bus.inputA_rad - fixtolerance
            candidateubound = bus.inputA_rad + fixtolerance
            lbound = max(lbound, candidatelbound)
            ubound = min(ubound, candidateubound)

        thetavar[bus] = model.addVar(
            obj=0.0, lb=lbound, ub=ubound, name="theta_" + str(bus.nodeID)
        )
        bus.thetavarind = varcount + newvarcount + 1
        newvarcount += 2

    cosvar = {}
    sinvar = {}
    thetaftvar = {}
    vfvtvar = {}

    logging.info("    Assumption. Phase angle diffs between -pi and pi.")

    for j in range(1, 1 + numbranches):
        branch = branches[j]
        f = branch.f
        t = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]
        busf = buses[count_of_f]
        bust = buses[count_of_t]

        cosvar[branch] = model.addVar(
            obj=0.0,
            lb=-1.0,
            ub=1.0,
            name="cos_%d_%d_%d" % (j, busf.nodeID, bust.nodeID),
        )
        sinvar[branch] = model.addVar(
            obj=0.0,
            lb=-1.0,
            ub=1.0,
            name="sin_%d_%d_%d" % (j, busf.nodeID, bust.nodeID),
        )

        # Major assumption! Heuristic: angle diffs between -pi and pi.
        thetaftvar[branch] = model.addVar(
            obj=0.0,
            lb=-math.pi,
            ub=math.pi,
            name="thetaft_%d_%d_%d" % (j, busf.nodeID, bust.nodeID),
        )

        ubound = busf.Vmax * bust.Vmax
        lbound = busf.Vmin * bust.Vmin

        vfvtvar[branch] = model.addVar(
            obj=0.0,
            lb=lbound,
            ub=ubound,
            name="vfvt_%d_%d_%d" % (j, busf.nodeID, bust.nodeID),
        )
        newvarcount += 4

    # Save polar variables data
    alldata["LP"]["vvar"] = vvar
    alldata["LP"]["thetavar"] = thetavar
    alldata["LP"]["cosvar"] = cosvar
    alldata["LP"]["sinvar"] = sinvar
    alldata["LP"]["thetaftvar"] = thetaftvar
    alldata["LP"]["vfvtvar"] = vfvtvar

    logging.info(
        "    Added %d new variables to handle polar formulation." % newvarcount
    )

    varcount += newvarcount


def lpformulator_ac_create_efvars(alldata, model, varcount):
    """Create nonconvex e, f variables for ACOPF"""

    evar = {}
    fvar = {}
    numbuses = alldata["numbuses"]
    buses = alldata["buses"]
    IDtoCountmap = alldata["IDtoCountmap"]
    efvarcount = 0
    fixtolerance = 1e-05
    if alldata["fixtolerance"] > 0:
        fixtolerance = alldata["fixtolerance"]

    logging.info("  Creating e,f variables.")

    for j in range(1, numbuses + 1):
        bus = buses[j]
        ubound = bus.Vmax
        lbound = -ubound

        if alldata["usemaxdispersion"]:
            lbound = bus.Vmin * math.cos(alldata["maxdispersion_rad"])
            ubound = bus.Vmax

        evar[bus] = model.addVar(
            obj=0.0, lb=lbound, ub=ubound, name="e_%d" % bus.nodeID
        )
        bus.evarind = varcount + efvarcount
        efvarcount += 1

        if alldata["usemaxdispersion"]:
            lbound = 0
            ubound = bus.Vmax * math.sin(alldata["maxdispersion_rad"])
        elif j == alldata["refbus"]:
            ubound = lbound = 0

        fvar[bus] = model.addVar(
            obj=0.0, lb=lbound, ub=ubound, name="f_%d" % bus.nodeID
        )
        bus.fvarind = varcount + efvarcount
        efvarcount += 1

    # Save e, f variables data
    alldata["LP"]["evar"] = evar
    alldata["LP"]["fvar"] = fvar

    logging.info("  Added %d e, f variables." % efvarcount)

    varcount += efvarcount


def lpformulator_ac_create_constraints(alldata, model):
    """Create constraints for ACOPF"""

    numbuses = alldata["numbuses"]
    buses = alldata["buses"]
    numbranches = alldata["numbranches"]
    branches = alldata["branches"]
    gens = alldata["gens"]
    IDtoCountmap = alldata["IDtoCountmap"]
    cvar = alldata["LP"]["cvar"]
    svar = alldata["LP"]["svar"]
    Pvar_f = alldata["LP"]["Pvar_f"]
    Pvar_t = alldata["LP"]["Pvar_t"]
    Qvar_f = alldata["LP"]["Qvar_f"]
    Qvar_t = alldata["LP"]["Qvar_t"]
    if alldata["branchswitching_mip"]:
        twinPvar_f = alldata["LP"]["twinPvar_f"]
        twinPvar_t = alldata["LP"]["twinPvar_t"]
        twinQvar_f = alldata["LP"]["twinQvar_f"]
        twinQvar_t = alldata["LP"]["twinQvar_t"]
    Pinjvar = alldata["LP"]["Pinjvar"]
    Qinjvar = alldata["LP"]["Qinjvar"]
    GenPvar = alldata["LP"]["GenPvar"]
    GenQvar = alldata["LP"]["GenQvar"]
    lincostvar = alldata["LP"]["lincostvar"]
    zvar = alldata["LP"]["zvar"]

    logging.info("Creating constraints.")
    logging.info("  Adding cost definition.")

    coeff = [gen.costvector[gen.costdegree - 1] for gen in gens.values()]
    variables = [GenPvar[gen] for gen in gens.values()]
    expr = gp.LinExpr(coeff, variables)
    model.addConstr(expr == lincostvar, name="lincostdef")

    numquadgens = 0
    for gen in gens.values():
        if gen.costdegree >= 2 and gen.costvector[0] > 0 and gen.status:
            numquadgens += 1

    logging.info(
        "    Number of generators with quadratic cost coefficient: %d." % numquadgens
    )

    if numquadgens > 0:
        if alldata["usequadcostvar"]:
            quadcostvar = alldata["LP"]["quadcostvar"]
            logging.info("    Adding quadcost definition constraint.")
            qcost = gp.QuadExpr()
            for gen in gens.values():
                if gen.costdegree == 2 and gen.costvector[0] != 0:
                    qcost.add(gen.costvector[0] * GenPvar[gen] * GenPvar[gen])

            model.addConstr(qcost <= quadcostvar, name="qcostdef")
        else:
            logging.info("    Adding quad cost to objective.")
            model.update()  # Necessary to flush changes in the objective function
            oldobj = model.getObjective()
            newobj = gp.QuadExpr(oldobj)
            for gen in gens.values():
                if gen.costdegree == 2 and gen.costvector[0] != 0:
                    newobj.add(gen.costvector[0] * GenPvar[gen] * GenPvar[gen])

            model.setObjective(newobj, GRB.MINIMIZE)

    # Active PF defs
    logging.info("  Adding active power flow definitions.")
    count = 0
    for j in range(1, 1 + numbranches):
        branch = branches[j]
        f = branch.f
        t = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]
        busf = buses[count_of_f]
        bust = buses[count_of_t]
        branch.Pfcname = "Pdef_%d_%d_%d" % (j, f, t)
        branch.Ptcname = "Pdef_%d_%d_%d" % (j, t, f)
        if branch.status:
            if alldata["substitute_nonconv"] == False or alldata["use_ef"] == False:
                # Gff cff + Gft cft + Bft sft
                expr = gp.LinExpr(
                    [branch.Gff, branch.Gft, branch.Bft],
                    [cvar[busf], cvar[branch], svar[branch]],
                )
                expP = Pvar_f[branch]
                if alldata["branchswitching_mip"]:
                    expP += twinPvar_f[branch]
                branch.Pdeffconstr = model.addConstr(expr == expP, name=branch.Pfcname)

                if alldata["branchswitching_mip"]:
                    if branch.constrainedflow:
                        coeff = branch.limit
                    else:
                        coeff = (
                            2 * alldata["sumPd"]
                        )  # the 2 is gratuitous but is an assumption nonetheless

                    model.addConstr(
                        Pvar_f[branch] <= coeff * zvar[branch],
                        name="upmip_P_%d_%d_%d" % (j, f, t),
                    )
                    model.addConstr(
                        Pvar_f[branch] >= -coeff * zvar[branch],
                        name="dnmip_P_%d_%d_%d" % (j, f, t),
                    )
                    model.addConstr(
                        twinPvar_f[branch] <= coeff * (1 - zvar[branch]),
                        name="upmip_twin_P_%d_%d_%d" % (j, f, t),
                    )
                    model.addConstr(
                        twinPvar_f[branch] >= -coeff * (1 - zvar[branch]),
                        name="dnmip_twin_P_%d_%d_%d" % (j, f, t),
                    )

                    count += 4

                # Gtt ctt + Gtf cft + Btf stf = Gtt ctt + Gtf cft - Btf sft
                expr = gp.LinExpr(
                    [
                        branch.Gtt,
                        branch.Gtf,
                        -branch.Btf,
                    ],  # minus because svarft = -svartf
                    [cvar[bust], cvar[branch], svar[branch]],
                )

                expP = Pvar_t[branch]
                if alldata["branchswitching_mip"]:
                    expP += twinPvar_t[branch]

                branch.Pdeftconstr = model.addConstr(expr == expP, name=branch.Ptcname)
                if alldata["branchswitching_mip"]:
                    if branch.constrainedflow:
                        coeff = branch.limit
                    else:
                        coeff = (
                            2 * alldata["sumPd"]
                        )  # the 2 is gratuitous but is an assumption nonetheless

                    model.addConstr(
                        Pvar_t[branch] <= coeff * zvar[branch],
                        name="upmip_P_%d_%d_%d" % (j, t, f),
                    )
                    model.addConstr(
                        Pvar_t[branch] >= -coeff * zvar[branch],
                        name="dnmip_P_%d_%d_%d" % (j, t, f),
                    )
                    model.addConstr(
                        twinPvar_t[branch] <= coeff * (1 - zvar[branch]),
                        name="upmip_twin_P_%d_%d_%d" % (j, t, f),
                    )
                    model.addConstr(
                        twinPvar_t[branch] >= -coeff * (1 - zvar[branch]),
                        name="dnmip_twin_P_%d_%d_%d" % (j, t, f),
                    )

                    count += 4

                count += 2
        else:
            branch.Pdeffconstr = model.addConstr(
                Pvar_f[branch] == 0, name=branch.Pfcname
            )
            branch.Pdeftconstr = model.addConstr(
                Pvar_t[branch] == 0, name=branch.Ptcname
            )

    logging.info("    %d active power flow definitions added." % count)

    # Reactive PF defs
    logging.info("  Adding reactive power flow definitions.")
    count = 0
    for j in range(1, 1 + numbranches):
        branch = branches[j]
        f = branch.f
        t = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]
        busf = buses[count_of_f]
        bust = buses[count_of_t]
        branch.Qfcname = "Qdef_%d_%d_%d" % (j, f, t)
        branch.Qtcname = "Qdef_%d_%d_%d" % (j, t, f)

        if branch.status:
            if alldata["substitute_nonconv"] == False or alldata["use_ef"] == False:
                # -Bff cff - Bft cft + Gft sft
                expr = gp.LinExpr(
                    [-branch.Bff, -branch.Bft, branch.Gft],
                    [cvar[busf], cvar[branch], svar[branch]],
                )

                expQ = Qvar_f[branch]
                if alldata["branchswitching_mip"]:
                    expQ += twinQvar_f[branch]
                branch.Qdeffconstr = model.addConstr(expr == expQ, name=branch.Qfcname)

                if alldata["branchswitching_mip"]:
                    if branch.constrainedflow:
                        coeff = branch.limit
                    else:
                        coeff = (
                            2 * alldata["sumQd"]
                        )  # the 2 is gratuitous but is an assumption nonetheless

                    model.addConstr(
                        Qvar_f[branch] <= coeff * zvar[branch],
                        name="upmip_Q_%d_%d_%d" % (j, f, t),
                    )
                    model.addConstr(
                        Qvar_f[branch] >= -coeff * zvar[branch],
                        name="dnmip_Q_%d_%d_%d" % (j, f, t),
                    )
                    model.addConstr(
                        twinQvar_f[branch] <= coeff * (1 - zvar[branch]),
                        name="upmip_twin_Q_%d_%d_%d" % (j, f, t),
                    )
                    model.addConstr(
                        twinQvar_f[branch] >= -coeff * (1 - zvar[branch]),
                        name="dnmip_twin_Q_%d_%d_%d" % (j, f, t),
                    )

                    count += 4

                # -Btt ctt - Btf cft + Gtf stf = -Btt ctt - Btf cft - Gtf sft
                expr = gp.LinExpr(
                    [-branch.Btt, -branch.Btf, -branch.Gtf],  # again, same minus
                    [cvar[bust], cvar[branch], svar[branch]],
                )

                expQ = Qvar_t[branch]
                if alldata["branchswitching_mip"]:
                    expQ += twinQvar_t[branch]
                branch.Qdeftconstr = model.addConstr(expr == expQ, name=branch.Qtcname)

                if alldata["branchswitching_mip"]:
                    if branch.constrainedflow:
                        coeff = branch.limit
                    else:
                        coeff = (
                            2 * alldata["sumQd"]
                        )  # the 2 is gratuitous but is an assumption nonetheless

                    model.addConstr(
                        Qvar_t[branch] <= coeff * zvar[branch],
                        name="upmip_Q_%d_%d_%d" % (j, t, f),
                    )
                    model.addConstr(
                        Qvar_t[branch] >= -coeff * zvar[branch],
                        name="dnmip_Q_%d_%d_%d" % (j, t, f),
                    )
                    model.addConstr(
                        twinQvar_t[branch] <= coeff * (1 - zvar[branch]),
                        name="upmip_twin_Q_%d_%d_%d" % (j, t, f),
                    )
                    model.addConstr(
                        twinQvar_t[branch] >= -coeff * (1 - zvar[branch]),
                        name="dnmip_twin_Q_%d_%d_%d" % (j, t, f),
                    )

                    count += 4

                count += 2
        else:  # out of operation
            branch.Qdeffconstr = model.addConstr(
                Qvar_f[branch] == 0, name=branch.Qfcname
            )
            branch.Qdeftconstr = model.addConstr(
                Qvar_t[branch] == 0, name=branch.Qtcname
            )

    logging.info("    %d reactive power flow definitions added." % count)

    # Balance constraints
    logging.info(
        "  Adding constraints stating bus injection = total outgoing power flow."
    )
    count = 0
    for j in range(1, 1 + numbuses):
        bus = buses[j]
        expr = gp.LinExpr()
        qexpr = gp.QuadExpr()
        doquad = False
        if bus.Gs != 0:
            expr.add(bus.Gs * cvar[bus])

        if alldata["branchswitching_comp"] == False:
            for branchid in bus.frombranchids.values():
                expr.add(Pvar_f[branches[branchid]])

            for branchid in bus.tobranchids.values():
                expr.add(Pvar_t[branches[branchid]])
        else:
            doquad = True
            # product of binary and powerflow
            for branchid in bus.frombranchids.values():
                qexpr.add((1 - zvar[branches[branchid]]) * Pvar_f[branches[branchid]])

            for branchid in bus.tobranchids.values():
                qexpr.add((1 - zvar[branches[branchid]]) * Pvar_t[branches[branchid]])

        if doquad:
            model.addConstr(
                expr + qexpr == Pinjvar[bus], name="PBaldef%d_%d" % (j, bus.nodeID)
            )
        else:
            model.addConstr(expr == Pinjvar[bus], name="PBaldef%d_%d" % (j, bus.nodeID))

        count += 1

    for j in range(1, 1 + numbuses):
        bus = buses[j]
        expr = gp.LinExpr()

        if bus.Bs != 0:
            expr.add(-bus.Bs * cvar[bus])
        if alldata["branchswitching_comp"] == False:
            for branchid in bus.frombranchids.values():
                expr.add(Qvar_f[branches[branchid]])

            for branchid in bus.tobranchids.values():
                expr.add(Qvar_t[branches[branchid]])

            model.addConstr(expr == Qinjvar[bus], name="QBaldef%d_%d" % (j, bus.nodeID))
        else:
            # product of binary and powerflow
            qexpr = gp.QuadExpr()

            for branchid in bus.frombranchids.values():
                qexpr.add((1 - zvar[branches[branchid]]) * Qvar_f[branches[branchid]])

            for branchid in bus.tobranchids.values():
                qexpr.add((1 - zvar[branches[branchid]]) * Qvar_t[branches[branchid]])

            model.addConstr(
                expr + qexpr == Qinjvar[bus], name="QBaldef%d_%d" % (j, bus.nodeID)
            )
        count += 1

    logging.info("    %d balance constraints added." % count)

    # Injection defs
    logging.info("  Adding injection definition constraints.")
    count = 0
    for j in range(1, 1 + numbuses):
        bus = buses[j]
        expr = gp.LinExpr()

        if len(bus.genidsbycount) > 0:
            for genid in bus.genidsbycount:
                gen = gens[genid]
                expr.add(GenPvar[gen])

        model.addConstr(Pinjvar[bus] == expr - bus.Pd, name="Bus_PInj_%d" % j)

        expr = gp.LinExpr()

        if len(bus.genidsbycount) > 0:
            for genid in bus.genidsbycount:
                gen = gens[genid]
                expr.add(GenQvar[gen])

        model.addConstr(Qinjvar[bus] == expr - bus.Qd, name="Bus_QInj_%d" % j)
        count += 2

    logging.info("    %d injection definition constraints added." % count)

    # Branch limits
    logging.info("  Adding branch limits.")
    count = 0
    for j in range(1, 1 + numbranches):
        branch = branches[j]

        if branch.status and branch.unboundedlimit == False:
            f = branch.f
            t = branch.t
            count_of_f = IDtoCountmap[f]
            count_of_t = IDtoCountmap[t]
            busf = buses[count_of_f]
            bust = buses[count_of_t]
            model.addConstr(
                Pvar_f[branch] * Pvar_f[branch] + Qvar_f[branch] * Qvar_f[branch]
                <= branch.limit**2,
                name="limit_f_%d_%d_%d" % (j, f, t),
            )
            model.addConstr(
                Pvar_t[branch] * Pvar_t[branch] + Qvar_t[branch] * Qvar_t[branch]
                <= branch.limit**2,
                name="limit_t_%d_%d_%d" % (j, t, f),
            )
            count += 2

    logging.info("    %d branch limits added." % count)

    # JABR.
    if alldata["skipjabr"] == False:
        logging.info("  Adding Jabr constraints.")
        count = 0
        for j in range(1, 1 + numbranches):
            branch = branches[j]
            if branch.status:
                f = branch.f
                t = branch.t
                count_of_f = IDtoCountmap[f]
                count_of_t = IDtoCountmap[t]
                busf = buses[count_of_f]
                bust = buses[count_of_t]
                model.addConstr(
                    cvar[branch] * cvar[branch] + svar[branch] * svar[branch]
                    <= cvar[busf] * cvar[bust],
                    name="jabr_%d_%d_%d" % (j, f, t),
                )
                count += 1

        logging.info("    %d Jabr constraints added." % count)
    else:
        logging.info("  Skipping Jabr inequalities.")

    # Active loss constraints.
    if alldata["useactivelossineqs"] == True:
        logging.info("  Adding active loss constraints.")
        count = 0
        for j in range(1, 1 + numbranches):
            branch = branches[j]
            if branch.status:
                f = branch.f
                t = branch.t
                count_of_f = IDtoCountmap[f]
                count_of_t = IDtoCountmap[t]
                busf = buses[count_of_f]
                bust = buses[count_of_t]

                model.addConstr(
                    cvar[busf] + cvar[bust] >= 2 * cvar[branch],
                    name="csc_%d_%d_%d" % (j, f, t),
                )

                count += 1

        if alldata["branchswitching_mip"]:
            # If not MIP, inequality implied by the csc ineq above.
            for j in range(1, 1 + numbranches):
                branch = branches[j]
                if branch.status and (branch.nongaining):  # (branch.isacline == True):

                    f = branch.f
                    t = branch.t
                    count_of_f = IDtoCountmap[f]
                    count_of_t = IDtoCountmap[t]
                    busf = buses[count_of_f]
                    bust = buses[count_of_t]
                    model.addConstr(
                        Pvar_f[branch] + Pvar_t[branch] >= 0,
                        name="aLa_%d_%d_%d" % (j, f, t),
                    )

                    count += 1

                    # Reactive, also.

                    if -branch.Bff > 0 or -branch.Btt > 0:
                        alpha = min(-branch.Bff, -branch.Btt)
                        beta = branch.Bft + branch.Btf
                        delta = 2 * alpha - beta
                        gamma = -branch.Gft + branch.Gtf

                        # print(alpha,beta, delta)
                        if delta >= 0:
                            model.addConstr(
                                Qvar_f[branch] + Qvar_t[branch] >= gamma * svar[branch],
                                name="rLa_%d_%d_%d" % (j, f, t),
                            )

                            count += 1

        logging.info("    %d active loss constraints added." % count)
    else:
        logging.info("  Skipping active loss inequalities.")

    # Polar representation
    if alldata["dopolar"]:
        lpformulator_ac_add_polarconstraints(alldata, model)

    # nonconvex e, f representation
    if alldata["use_ef"] and alldata["useconvexformulation"] == False:
        lpformulator_ac_add_nonconvexconstraints(alldata, model)

    if alldata["branchswitching_mip"] or alldata["branchswitching_comp"]:
        boundzs = True
        if boundzs:
            exp = gp.LinExpr()
            delta = numbranches
            N = numbranches - delta  # <<<<<<---- here is the heuristic lower bound
            logging.info("In bound_zs constraint, N = %d." % (N))
            for j in range(1, 1 + numbranches):
                branch = branches[j]
                exp += zvar[branch]
            model.addConstr(exp >= N, name="sumz_lower")


def lpformulator_ac_add_polarconstraints(alldata, model):
    """Create polar representation constraints for ACOPF"""

    buses = alldata["buses"]
    numbuses = alldata["numbuses"]
    branches = alldata["branches"]
    numbranches = alldata["numbranches"]
    vvar = alldata["LP"]["vvar"]
    thetavar = alldata["LP"]["thetavar"]
    thetaftvar = alldata["LP"]["thetaftvar"]
    vfvtvar = alldata["LP"]["vfvtvar"]
    cosvar = alldata["LP"]["cosvar"]
    sinvar = alldata["LP"]["sinvar"]
    IDtoCountmap = alldata["IDtoCountmap"]
    cvar = alldata["LP"]["cvar"]
    svar = alldata["LP"]["svar"]
    count = 0

    logging.info("  Adding polar constraints.")

    for j in range(1, 1 + numbuses):
        bus = buses[j]
        model.addConstr(cvar[bus] == vvar[bus] ** 2, name="cffdef_%d" % j)
        count += 1

    for j in range(1, 1 + numbranches):
        branch = branches[j]
        f = branch.f
        t = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]
        busf = buses[count_of_f]
        bust = buses[count_of_t]
        model.addConstr(
            thetaftvar[branch] == thetavar[busf] - thetavar[bust],
            name="thetaftdef_%d" % j,
        )

        gc = model.addGenConstrCos(
            thetaftvar[branch], cosvar[branch], name="cosdef_%d_%d_%d" % (j, f, t)
        )
        gs = model.addGenConstrSin(
            thetaftvar[branch], sinvar[branch], name="sindef_%d_%d_%d" % (j, f, t)
        )
        model.addConstr(
            vfvtvar[branch] == vvar[busf] * vvar[bust],
            name="vfvtdef_%d_%d_%d" % (j, f, t),
        )

        model.addConstr(
            cvar[branch] == vfvtvar[branch] * cosvar[branch], name="cftdef_%d" % j
        )
        model.addConstr(
            svar[branch] == vfvtvar[branch] * sinvar[branch], name="sftdef_%d" % j
        )
        count += 6  # Count general constraints as well

    logging.info("    %d polar constraints added." % count)


def lpformulator_ac_add_nonconvexconstraints(alldata, model):
    """Create nonconvex e, f constraints"""

    buses = alldata["buses"]
    numbuses = alldata["numbuses"]
    branches = alldata["branches"]
    numbranches = alldata["numbranches"]
    evar = alldata["LP"]["evar"]
    fvar = alldata["LP"]["fvar"]
    cvar = alldata["LP"]["cvar"]
    svar = alldata["LP"]["svar"]
    IDtoCountmap = alldata["IDtoCountmap"]
    count = 0
    logging.info("  Adding nonconvex e, f, constraints.")

    for j in range(1, 1 + numbuses):
        bus = buses[j]
        model.addConstr(
            -cvar[bus] + evar[bus] * evar[bus] + fvar[bus] * fvar[bus] == 0,
            name="cbusdef_%d_%d" % (j, bus.nodeID),
        )
        count += 1

    for j in range(1, 1 + numbranches):
        branch = branches[j]

        if branch.status:
            f = branch.f
            t = branch.t
            count_of_f = IDtoCountmap[f]
            count_of_t = IDtoCountmap[t]
            busf = buses[count_of_f]
            bust = buses[count_of_t]

            model.addConstr(
                -cvar[branch] + evar[busf] * evar[bust] + fvar[busf] * fvar[bust] == 0,
                name="cdef_%d_%d_%d" % (j, f, t),
            )
            model.addConstr(
                -svar[branch] - evar[busf] * fvar[bust] + fvar[busf] * evar[bust] == 0,
                name="sdef_%d_%d_%d" % (j, f, t),
            )
            count += 2

    logging.info("    %d nonconvex e, f constraints added." % count)


def grbderive_xtra_sol_values_fromvoltages(alldata):
    # Generates complete solution vectors from input voltages.

    model = alldata["model"]

    model.update()

    numvars = model.NumVars
    xbuffer = alldata["LP"]["xbuffer"] = {}  # dictionary to store solution values

    numbuses = alldata["numbuses"]
    buses = alldata["buses"]
    IDtoCountmap = alldata["IDtoCountmap"]

    branches = alldata["branches"]
    numbranches = alldata["numbranches"]

    cvar = alldata["LP"]["cvar"]
    svar = alldata["LP"]["svar"]

    logging.info("Creating derived solution values from input voltage solution.")

    loud = False

    for j in range(1, numbuses + 1):
        bus = buses[j]
        bus.inpute = bus.inputV * math.cos(bus.inputA_rad)
        bus.inputf = bus.inputV * math.sin(bus.inputA_rad)

    if alldata["use_ef"]:
        for j in range(1, numbuses + 1):
            bus = buses[j]
            xbuffer[evar[bus]] = bus.inpute
            xbuffer[fvar[bus]] = bus.inputf

    if alldata["dopolar"]:
        vvar = alldata["LP"]["vvar"]
        thetavar = alldata["LP"]["vvar"]
        for j in range(1, numbuses + 1):
            bus = buses[j]
            xbuffer[vvar[bus]] = bus.inputV
            xbuffer[
                thetavar[bus]
            ] = bus.inputA_rad  # here is the bug; the first 'bus.' should not be there

    if alldata["dopolar"] == False:

        for j in range(1, numbuses + 1):
            bus = buses[j]
            xbuffer[cvar[bus]] = bus.inpute * bus.inpute + bus.inputf * bus.inputf

        for j in range(1, 1 + numbranches):
            branch = branches[j]
            if branch.status:
                f = branch.f
                t = branch.t
                count_of_f = IDtoCountmap[f]
                count_of_t = IDtoCountmap[t]
                busf = buses[count_of_f]
                bust = buses[count_of_t]

                xbuffer[cvar[branch]] = (
                    busf.inpute * bust.inpute + busf.inputf * bust.inputf
                )
                xbuffer[svar[branch]] = (
                    -busf.inpute * bust.inputf + busf.inputf * bust.inpute
                )

    else:
        # Note: we may not need all of thse
        vfvtvar = alldata["LP"]["vfvtvar"]
        thetaftvar = alldata["LP"]["thetaftvar"]
        cosvar = alldata["LP"]["cosvar"]
        sinvar = alldata["LP"]["sinvar"]

        for j in range(1, numbuses + 1):
            bus = buses[j]
            xbuffer[cvar[bus]] = bus.inputV * bus.inputV

        for j in range(1, 1 + numbranches):
            branch = branches[j]
            if branch.status:
                f = branch.f
                t = branch.t
                count_of_f = IDtoCountmap[f]
                count_of_t = IDtoCountmap[t]
                busf = buses[count_of_f]
                bust = buses[count_of_t]
                vfvt = xbuffer[vfvtvar[branch]] = busf.inputV * bust.inputV
                thetaft = xbuffer[thetaftvar[branch]] = (
                    busf.inputA_rad - bust.inputA_rad
                )
                c = xbuffer[cosvar[branch]] = math.cos(thetaft)
                s = xbuffer[sinvar[branch]] = math.sin(thetaft)
                xbuffer[cvar[branch]] = vfvt * c
                xbuffer[svar[branch]] = vfvt * s

    logging.info("Derived e, f, c, s values.")

    Pvar_f = alldata["LP"]["Pvar_f"]
    Pvar_t = alldata["LP"]["Pvar_t"]
    Qvar_f = alldata["LP"]["Qvar_f"]
    Qvar_t = alldata["LP"]["Qvar_t"]

    for j in range(1, 1 + numbranches):
        branch = branches[j]
        row = model.getRow(branch.Pdeffconstr)
        if (
            loud
        ):  # TODO-Dan We have to discuss the loud variable. If you want to keep it then it should be handled
            # properly via a setting and not be hard-coded
            print(
                j,
                branch.Pdeffconstr,
                row,
                branch.Pdeffconstr.Sense,
                branch.Pdeffconstr.RHS,
            )
            print(row.size())
        sum = -branch.Pdeffconstr.RHS
        leadcoeff = 0
        for i in range(row.size()):
            var = row.getVar(i)
            coeff = row.getCoeff(i)
            if var.Varname != Pvar_f[branch].Varname:
                if loud:
                    print("   ", i, coeff, var.Varname, "at", xbuffer[var])
                sum += coeff * xbuffer[var]
            else:
                if loud:
                    print("   >", i, coeff, var.Varname)
                leadcoeff = coeff
        if loud:
            print("sum =", sum, leadcoeff)
        xbuffer[Pvar_f[branch]] = -sum / leadcoeff
        # leadcoeff should be +1 or -1
        logging.info(
            "%s derived at %f." % (Pvar_f[branch].Varname, xbuffer[Pvar_f[branch]])
        )

        row = model.getRow(branch.Pdeftconstr)
        if loud:
            print(
                j,
                branch.Pdeftconstr,
                row,
                branch.Pdeftconstr.Sense,
                branch.Pdeftconstr.RHS,
            )
            print(row.size())
        sum = -branch.Pdeftconstr.RHS
        leadcoeff = 0
        for i in range(row.size()):
            var = row.getVar(i)
            coeff = row.getCoeff(i)
            if var.Varname != Pvar_t[branch].Varname:
                if loud:
                    print("   ", i, coeff, var.Varname, "at", xbuffer[var])
                sum += coeff * xbuffer[var]
            else:
                if loud:
                    print("   >", i, coeff, var.Varname)
                leadcoeff = coeff
        if loud:
            print("sum =", sum, leadcoeff)
        xbuffer[Pvar_t[branch]] = -sum / leadcoeff
        # leadcoeff should be +1 or -1
        logging.info(
            "%s derived at %f." % (Pvar_t[branch].Varname, xbuffer[Pvar_t[branch]])
        )

        ##################### Now, reactive power flows

        row = model.getRow(branch.Qdeffconstr)
        if loud:
            print(
                j,
                branch.Qdeffconstr,
                row,
                branch.Qdeffconstr.Sense,
                branch.Qdeffconstr.RHS,
            )
            print("size", row.size())
        sum = -branch.Qdeffconstr.RHS
        leadcoeff = 0

        for i in range(row.size()):
            var = row.getVar(i)
            coeff = row.getCoeff(i)
            if var.Varname != Qvar_f[branch].Varname:
                if loud:
                    print("   ", i, coeff, var.Varname, "at", xbuffer[var])
                sum += coeff * xbuffer[var]
            else:
                if loud:
                    print("   >", i, coeff, var.Varname)
                leadcoeff = coeff
        if loud:
            print("sum =", sum, leadcoeff)

        xbuffer[Qvar_f[branch]] = -sum / leadcoeff
        # leadcoeff should be +1 or -1
        logging.info(
            "%s derived at %f." % (Qvar_f[branch].Varname, xbuffer[Qvar_f[branch]])
        )

        row = model.getRow(branch.Qdeftconstr)
        if loud:
            print(
                "\n",
                j,
                branch.Qdeftconstr,
                row,
                branch.Qdeftconstr.Sense,
                branch.Qdeftconstr.RHS,
            )
            print(row.size())
        sum = -branch.Qdeftconstr.RHS
        leadcoeff = 0
        for i in range(row.size()):
            var = row.getVar(i)
            coeff = row.getCoeff(i)
            if var.Varname != Qvar_t[branch].Varname:
                product = coeff * xbuffer[var]
                sum += product
                if loud:
                    print(
                        "   ",
                        i,
                        coeff,
                        var.Varname,
                        "at",
                        xbuffer[var],
                        " prod:",
                        product,
                    )
                    print("    sum:", sum)

            else:
                if loud:
                    print("   >", i, coeff, var.Varname)
                leadcoeff = coeff
        if loud:
            print("sum =", sum, leadcoeff)
        xbuffer[Qvar_t[branch]] = -sum / leadcoeff
        # leadcoeff should be +1 or -1
        logging.info(
            "%s derived at %f." % (Qvar_t[branch].Varname, xbuffer[Qvar_t[branch]])
        )

    # next, power flow injections
    Pinjvar = alldata["LP"]["Pinjvar"]
    for j in range(1, 1 + numbuses):
        bus = buses[j]

        injection = 0
        for branchid in bus.frombranchids.values():
            branch = branches[branchid]
            injection += xbuffer[Pvar_f[branch]]

        for branchid in bus.tobranchids.values():
            branch = branches[branchid]
            injection += xbuffer[Pvar_t[branch]]

        xbuffer[Pinjvar[bus]] = injection

    # next, power flow injections
    Qinjvar = alldata["LP"]["Qinjvar"]
    for j in range(1, 1 + numbuses):
        bus = buses[j]

        injection = 0
        for branchid in bus.frombranchids.values():
            branch = branches[branchid]
            injection += xbuffer[Qvar_f[branch]]

        for branchid in bus.tobranchids.values():
            branch = branches[branchid]
            injection += xbuffer[Qvar_t[branch]]

        xbuffer[Qinjvar[bus]] = injection


def lpformulator_ac_strictchecker(alldata, model, spitoutvector):
    """Checks feasibility of input solution -- reports infeasibilities"""
    logging.info("Strict feasibility check for input voltage solution.\n")

    if alldata["strictcheckvoltagesolution"]:
        grbderive_xtra_sol_values_fromvoltages(alldata)

    buses = alldata["buses"]
    numbuses = alldata["numbuses"]
    branches = alldata["branches"]
    numbranches = alldata["numbranches"]
    gens = alldata["gens"]

    max_violation_string = None
    max_violation_value = 0

    alldata["violation"] = {}  # dictionary to keep track of violations
    Vmagviol = alldata["violation"]["Vmagviol"] = {}
    IPviol = alldata["violation"]["IPviol"] = {}
    IQviol = alldata["violation"]["IQviol"] = {}
    branchlimitviol = alldata["violation"]["branchlimit"] = {}
    for j in range(1, numbuses + 1):
        bus = buses[j]
        alldata["violation"][bus] = {}

    if alldata["use_ef"]:
        evar = alldata["LP"]["evar"]
        fvar = alldata["LP"]["fvar"]

    xbuffer = alldata["LP"]["xbuffer"]

    maxlbviol = maxubviol = 0
    badlbvar = None
    badubvar = None

    logging.info("Direct bus magnitude bounds check. Error issued if infeasibility.")

    for j in range(1, numbuses + 1):
        bus = buses[j]
        if bus.inputV > bus.Vmax:
            logging.error(
                ">>> Error: bus # %d has input voltage %f which is larger than Vmax %f."
                % (j, bus.inputV, bus.Vmax)
            )
            thisviol = bus.inputV - bus.Vmax
            if thisviol > max_violation_value:
                max_violation_string = "bus_" + str(bus.nodeID) + "_Vmax"
                max_violation_value = thisviol
        if bus.inputV < bus.Vmin:
            logging.error(
                ">>> Error: bus # %d has input voltage %f which is smaller than Vmin %f."
                % (j, bus.inputV, bus.Vmin)
            )
            thisviol = bus.Vmin - bus.inputV
            if thisviol > max_violation_value:
                max_violation_string = "bus_" + str(bus.nodeID) + "_Vmin"
                max_violation_value = thisviol

        alldata["violation"][bus]["Vmax"] = max(bus.inputV - bus.Vmax, 0)
        alldata["violation"][bus]["Vmin"] = max(bus.Vmin - bus.inputV, 0)

        candmaxviol = alldata["violation"][bus]["Vmax"]
        if candmaxviol < alldata["violation"][bus]["Vmin"]:
            candmaxviol = -alldata["violation"][bus]["Vmin"]

        Vmagviol[bus] = candmaxviol

    logging.info("Checked input Vmag values.")

    logging.info("Direct branch limit check. Error issued if infeasible.")

    Pvar_f = alldata["LP"]["Pvar_f"]
    Pvar_t = alldata["LP"]["Pvar_t"]
    Qvar_f = alldata["LP"]["Qvar_f"]
    Qvar_t = alldata["LP"]["Qvar_t"]

    for j in range(1, 1 + numbranches):
        branch = branches[j]
        fromvalue = math.sqrt(
            xbuffer[Pvar_f[branch]] * xbuffer[Pvar_f[branch]]
            + xbuffer[Qvar_f[branch]] * xbuffer[Qvar_f[branch]]
        )
        fromviol = max(fromvalue - branch.limit, 0)
        if fromvalue > branch.limit:
            logging.error(
                ">>> Error: branch # %d has 'from' flow magnitude %f which is larger than limit %f."
                % (j, fromvalue, branch.limit)
            )
            logging.error("    branch is ( %d %d )." % (branch.f, branch.t))
            thisviol = fromviol
            if thisviol > max_violation_value:
                max_violation_string = "branch_" + str(j) + "_from"
                max_violation_value = thisviol

        tovalue = math.sqrt(
            xbuffer[Pvar_t[branch]] * xbuffer[Pvar_t[branch]]
            + xbuffer[Qvar_t[branch]] * xbuffer[Qvar_t[branch]]
        )
        toviol = max(tovalue - branch.limit, 0)
        if tovalue > branch.limit:
            logging.error(
                ">>> Error: branch # %d has 'to' flow magnitude %f which is larger than limit %f."
                % (j, tovalue, branch.limit)
            )
            logging.error("    branch is ( %d %d )." % (branch.f, branch.t))
            thisviol = toviol
            if thisviol > max_violation_value:
                max_violation_string = "branch_" + str(j) + "_to"
                max_violation_value = thisviol
        alldata["violation"]["branchlimit"][branch] = max(fromviol, toviol)

    if alldata["use_ef"]:
        for j in range(1, numbuses + 1):
            bus = buses[j]
            (
                maxlbviol,
                maxubviol,
                badlbvar,
                badubvar,
                max_violation_value,
                max_violation_string,
            ) = lpformulator_checkviol_simple(
                alldata,
                model,
                evar[bus],
                bus.inpute,
                maxlbviol,
                maxubviol,
                badlbvar,
                badubvar,
                max_violation_value,
                max_violation_string,
                True,
            )

            (
                maxlbviol,
                maxubviol,
                badlbvar,
                badubvar,
                max_violation_value,
                max_violation_string,
            ) = lpformulator_checkviol_simple(
                alldata,
                model,
                fvar[bus],
                bus.inputf,
                maxlbviol,
                maxubviol,
                badlbvar,
                badubvar,
                max_violation_value,
                max_violation_string,
                True,
            )

            alldata["LP"]["xbuffer"][evar[bus]] = bus.inpute
            alldata["LP"]["xbuffer"][fvar[bus]] = bus.inputf

        logging.info("e, f values.")

    if alldata["dopolar"]:
        logging.info(
            "Checking polar quantities. Note: bounds shifted by input solution."
        )  # which will repeat the Vmag value check, so ...
        vvar = alldata["LP"]["vvar"]
        thetavar = alldata["LP"]["vvar"]
        for j in range(1, numbuses + 1):
            bus = buses[j]
            (
                maxlbviol,
                maxubviol,
                badlbvar,
                badubvar,
                max_violation_value,
                max_violation_string,
            ) = lpformulator_checkviol_simple(
                alldata,
                model,
                vvar[bus],
                bus.inputV,
                maxlbviol,
                maxubviol,
                badlbvar,
                badubvar,
                max_violation_value,
                max_violation_string,
                True,
            )

        logging.info("Vmag values checked.")

    logging.info("Checking power flow values.")

    for j in range(1, 1 + numbranches):
        branch = branches[j]
        (
            maxlbviol,
            maxubviol,
            badlbvar,
            badubvar,
            max_violation_value,
            max_violation_string,
        ) = lpformulator_checkviol_simple(
            alldata,
            model,
            Pvar_f[branch],
            xbuffer[Pvar_f[branch]],
            maxlbviol,
            maxubviol,
            badlbvar,
            badubvar,
            max_violation_value,
            max_violation_string,
            True,
        )
        (
            maxlbviol,
            maxubviol,
            badlbvar,
            badubvar,
            max_violation_value,
            max_violation_string,
        ) = lpformulator_checkviol_simple(
            alldata,
            model,
            Pvar_t[branch],
            xbuffer[Pvar_t[branch]],
            maxlbviol,
            maxubviol,
            badlbvar,
            badubvar,
            max_violation_value,
            max_violation_string,
            True,
        )

    for j in range(1, 1 + numbranches):
        branch = branches[j]
        (
            maxlbviol,
            maxubviol,
            badlbvar,
            badubvar,
            max_violation_value,
            max_violation_string,
        ) = lpformulator_checkviol_simple(
            alldata,
            model,
            Qvar_f[branch],
            xbuffer[Qvar_f[branch]],
            maxlbviol,
            maxubviol,
            badlbvar,
            badubvar,
            max_violation_value,
            max_violation_string,
            True,
        )
        (
            maxlbviol,
            maxubviol,
            badlbvar,
            badubvar,
            max_violation_value,
            max_violation_string,
        ) = lpformulator_checkviol_simple(
            alldata,
            model,
            Qvar_t[branch],
            xbuffer[Qvar_t[branch]],
            maxlbviol,
            maxubviol,
            badlbvar,
            badubvar,
            max_violation_value,
            max_violation_string,
            True,
        )

    logging.info("Checking flow balance constraints.")
    Pinjvar = alldata["LP"]["Pinjvar"]
    Qinjvar = alldata["LP"]["Qinjvar"]

    for j in range(1, 1 + numbuses):
        bus = buses[j]

        (
            maxlbviol,
            maxubviol,
            badlbvar,
            badubvar,
            max_violation_value,
            max_violation_string,
        ) = lpformulator_checkviol_simple(
            alldata,
            model,
            Pinjvar[bus],
            xbuffer[Pinjvar[bus]],
            maxlbviol,
            maxubviol,
            badlbvar,
            badubvar,
            max_violation_value,
            max_violation_string,
            True,
        )

        alldata["violation"][bus]["Pinjmax"] = max(
            xbuffer[Pinjvar[bus]] - Pinjvar[bus].ub, 0
        )
        alldata["violation"][bus]["Pinjmin"] = max(
            Pinjvar[bus].lb - xbuffer[Pinjvar[bus]], 0
        )

        injection = 0
        for branchid in bus.frombranchids.values():
            branch = branches[branchid]
            injection += xbuffer[Pvar_f[branch]]

        for branchid in bus.tobranchids.values():
            branch = branches[branchid]
            injection += xbuffer[Pvar_t[branch]]

        # Pinjvar variable bounds should accommodate computed injection
        # get bounds for P and Q injection variables
        Pubound, Plbound, Qubound, Qlbound = computebalbounds(alldata, bus)

        # but also, construct min/max injection at the bus by looking at available generators and load
        myPubound = myPlbound = 0
        for gencounter in bus.genidsbycount:
            if gens[gencounter].status:
                myPubound += gens[gencounter].Pmax
                myPlbound += gens[gencounter].Pmin

        minnetgen = myPlbound - bus.Pd
        maxnetgen = myPubound - bus.Pd

        logging.info(
            "Bus ID %s #%d injection %g mingen %g maxgen %g load %g."
            % (bus.nodeID, j, injection, myPlbound, myPubound, bus.Pd)
        )
        logging.info("   min net generation %g; max %g." % (minnetgen, maxnetgen))

        candmaxviol = alldata["violation"][bus]["Pinjmax"]
        if candmaxviol < alldata["violation"][bus]["Pinjmin"]:
            candmaxviol = -alldata["violation"][bus]["Pinjmin"]
        IPviol[bus] = candmaxviol

    for j in range(1, 1 + numbuses):
        bus = buses[j]

        (
            maxlbviol,
            maxubviol,
            badlbvar,
            badubvar,
            max_violation_value,
            max_violation_string,
        ) = lpformulator_checkviol_simple(
            alldata,
            model,
            Qinjvar[bus],
            xbuffer[Qinjvar[bus]],
            maxlbviol,
            maxubviol,
            badlbvar,
            badubvar,
            max_violation_value,
            max_violation_string,
            True,
        )

        alldata["violation"][bus]["Qinjmax"] = max(
            xbuffer[Qinjvar[bus]] - Qinjvar[bus].ub, 0
        )
        alldata["violation"][bus]["Qinjmin"] = max(
            Qinjvar[bus].lb - xbuffer[Qinjvar[bus]], 0
        )

        injection = 0
        for branchid in bus.frombranchids.values():
            branch = branches[branchid]
            injection += xbuffer[Qvar_f[branch]]

        for branchid in bus.tobranchids.values():
            branch = branches[branchid]
            injection += xbuffer[Qvar_t[branch]]

        # Qinjvar variable bounds should accommodate computed injection
        # get bounds for P and Q injection variables
        Pubound, Plbound, Qubound, Qlbound = computebalbounds(alldata, bus)
        # print(Plbound, Pubound)

        # but also, construct min/max injection at the bus by looking at available generators and load
        myQubound = myQlbound = 0
        for gencounter in bus.genidsbycount:
            if gens[gencounter].status:
                myQubound += gens[gencounter].Qmax
                myQlbound += gens[gencounter].Qmin

        minnetgen = myQlbound - bus.Qd
        maxnetgen = myQubound - bus.Qd

        logging.info(
            "Bus ID %s #%d injection %g mingen %g maxgen %g load %g."
            % (bus.nodeID, j, injection, myQlbound, myQubound, bus.Qd)
        )
        logging.info("   min net generation %g; max %g." % (minnetgen, maxnetgen))

        candmaxviol = alldata["violation"][bus]["Qinjmax"]
        if candmaxviol < alldata["violation"][bus]["Qinjmin"]:
            candmaxviol = -alldata["violation"][bus]["Qinjmin"]
        IQviol[bus] = candmaxviol

    worstboundviol_report(badlbvar, maxlbviol, "LB")
    worstboundviol_report(badubvar, maxubviol, "UB")

    logging.info("\nSummary: Max LB viol %g, Max UB viol %g." % (maxlbviol, maxubviol))
    logging.info(
        "\nSummary: Max violation %g, key: %s."
        % (max_violation_value, max_violation_string)
    )

    print(alldata["dographics"])
    if alldata["dographics"]:
        textlist = []
        grbgraphical(alldata, "violation", textlist)

    break_exit("strictgraphical")


def lpformulator_checkviol_simple(
    alldata,
    model,
    grbvariable,
    value,
    maxlbviol,
    maxubviol,
    badlbvar,
    badubvar,
    max_violation_value,
    max_violation_string,
    loud,
):
    """
    Returns bounds infeasibility if setting grbvariable to value value
    maxlbviol, maxubviol, badlb, ubvar are updated if the infeasibility is larger
    """

    ub = grbvariable.ub
    lb = grbvariable.lb

    if loud:
        logging.info(
            "%s =  %.16e  [ LB %.4e  UB %.4e ]" % (grbvariable.varname, value, lb, ub)
        )

    lbviol = lb - value
    ubviol = value - ub

    if lbviol > 0:
        logging.info(
            "LBVIOL %s  LB %.6e  x %.16e  UB %.6e"
            % (grbvariable.varname, lb, value, ub)
        )

    if lbviol > maxlbviol:
        logging.info(" -> incumbent MAX LBVIOL")
        maxlbviol = lb - value
        badlbvar = grbvariable

    if ubviol > 0:
        logging.info(
            "UBVIOL %s  LB %.6e  x %.16e  UB %.6e"
            % (grbvariable.varname, lb, value, ub)
        )

    if ubviol > maxubviol:
        logging.info(" -> incumbent MAX UBVIOL")
        maxubviol = value - ub
        badubvar = grbvariable

    if maxubviol > max_violation_value:
        max_violation_string = grbvariable.Varname + "_ub"
        max_violation_value = maxubviol
    if maxlbviol > max_violation_value:
        max_violation_string = grbvariable.Varname + "_lb"
        max_violation_value = maxlbviol

    return (
        maxlbviol,
        maxubviol,
        badlbvar,
        badubvar,
        max_violation_value,
        max_violation_string,
    )


def worstboundviol_report(badvar, maxviol, boundtype):
    if badvar != None:
        logging.info(
            "Worst %s bound violation by %s viol %g."
            % (boundtype, badvar.Varname, maxviol)
        )
    else:
        logging.info("No %s bound violations.\n" % boundtype)


def writempsfile(alldata, model, filename):
    logging.info("Writing mpsfile to %s." % (filename))
    model.write(filename)
    break_exit("wrotefile")


def writemipstart(alldata):
    # Write z variables

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

    break_exit("wrote mipstart")
