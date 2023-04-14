import math
import time
import logging
import numpy as np
import gurobipy as gp
from gurobipy import GRB

from .utils import break_exit
from .grbfile import grbreadvoltsfile


def lpformulator_iv_body(alldata, model):
    """
    Adds variables and constraints for IV formulation to a given Gurobi model

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    model : gurobipy.Model
        Gurobi model to be constructed
    """

    logger = logging.getLogger("OpfLogger")
    # Create model variables
    lpformulator_iv_create_vars(alldata, model)
    # Create model constraints
    lpformulator_iv_create_constraints(alldata, model)

    model.update()  # Update to get correct model stats
    logger.info(
        "Constructed ACOPF model with %d variables and %d constraints.\n"
        % (model.NumVars, model.NumConstrs)
    )

    if alldata["lpfilename"] != None:
        model.write(
            alldata["lpfilename"]
        )  # FIXME remove.  Jarek: I am using this for debugging, for now
        logger.info("Wrote LP to " + alldata["lpfilename"])

    alldata["model"] = model


def lpformulator_iv_create_vars(alldata, model):
    """
    Creates and adds variables for IV formulation to a given Gurobi model

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    model : gurobipy.Model
        Gurobi model to be constructed
    """

    logger = logging.getLogger("OpfLogger")
    logger.info("Creating variables.")

    fixtolerance = 1e-05

    if alldata["fixtolerance"] > 0:
        fixtolerance = alldata["fixtolerance"]

    numbuses = alldata["numbuses"]
    buses = alldata["buses"]
    IDtoCountmap = alldata["IDtoCountmap"]
    gens = alldata["gens"]

    # cffvars
    cffvar = {}

    varcount = 0

    """
    for j in range(1,numbuses+1):
        bus = buses[j]
        # First, injection variables
        maxprod = bus.Vmax*bus.Vmax
        minprod = bus.Vmin*bus.Vmin
        ubound  = maxprod
        lbound  = minprod

        if alldata['fixcs'] and bus.inputvoltage:
            lbound = bus.inputV*bus.inputV - fixtolerance
            ubound = bus.inputV*bus.inputV + fixtolerance

        cffvar[bus] = model.addVar(obj = 0.0, lb = lbound, ub = ubound,
                                 name = "cff_%d"%(bus.nodeID))

        bus.cffvarind = varcount
        varcount += 1

        #at this point, minprod is the square of bus min voltage

    logger.info('Added %d cff variables\n'%(varcount))
    """

    lpformulator_iv_create_efvars(alldata, model, varcount)

    # Next, generator variables
    GenPvar = {}
    GenQvar = {}

    for j in range(1, numbuses + 1):
        bus = buses[j]
        # First, injection variables
        maxprod = bus.Vmax * bus.Vmax
        minprod = bus.Vmin * bus.Vmin

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

    logger.info("Added generator variables.")

    # Branch related variables
    branches = alldata["branches"]
    numbranches = alldata["numbranches"]

    irvar_f = {}  # branches!
    irvar_t = {}
    ijvar_f = {}  # branches!
    ijvar_t = {}

    Pvar_f = {}
    Qvar_f = {}
    Pvar_t = {}
    Qvar_t = {}

    if alldata["ivtype"] == "plain":
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
            if branch.constrainedflow:
                limit = branch.limit
            else:
                limit = 2 * (
                    abs(alldata["summaxgenP"]) + abs(alldata["summaxgenQ"])
                )  # Generous: assumes line charging up to 100%. However it still amounts to an assumption.

            ubound = 1e5  # 1 + np.sqrt(limit*limit/minprod ) #upper bound on current magnitude
            lbound = -ubound

            irvar_f[branch] = model.addVar(
                obj=0.0,
                lb=lbound,
                ub=ubound,
                name="ir_f_%d_%d_%d" % (j, busf.nodeID, bust.nodeID),
            )
            varcount += 1
            irvar_t[branch] = model.addVar(
                obj=0.0,
                lb=lbound,
                ub=ubound,
                name="ir_t_%d_%d_%d" % (j, bust.nodeID, busf.nodeID),
            )
            varcount += 1

            ijvar_f[branch] = model.addVar(
                obj=0.0,
                lb=lbound,
                ub=ubound,
                name="ij_f_%d_%d_%d" % (j, busf.nodeID, bust.nodeID),
            )
            varcount += 1
            ijvar_t[branch] = model.addVar(
                obj=0.0,
                lb=lbound,
                ub=ubound,
                name="ij_t_%d_%d_%d" % (j, bust.nodeID, busf.nodeID),
            )
            varcount += 1

        logger.info("Added branch current variables.")
    else:
        logger.info("Aggressive IV formulation, so skipping current variables.")

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

        Pvar_t[branch] = model.addVar(
            obj=0.0,
            lb=lbound,
            ub=ubound,
            name="P_%d_%d_%d" % (j, bust.nodeID, busf.nodeID),
        )
        branch.Ptfvarind = varcount
        varcount += 1

        Qvar_f[branch] = model.addVar(
            obj=0.0,
            lb=lbound,
            ub=ubound,
            name="Q_%d_%d_%d" % (j, busf.nodeID, bust.nodeID),
        )
        branch.Qftvarind = varcount
        varcount += 1

        Qvar_t[branch] = model.addVar(
            obj=0.0,
            lb=lbound,
            ub=ubound,
            name="Q_%d_%d_%d" % (j, bust.nodeID, busf.nodeID),
        )
        branch.Qtfvarind = varcount
        varcount += 1

    logger.info("Added branch power flow variables.")

    """
    #next, bus current variables
    #coarse bounds for now

    #commented out but will probably dump
    for j in range(1,numbuses+1):
        bus = buses[j]
        # First, injection variables
        ubound = 2*(abs(alldata['summaxgenP']) + abs(alldata['summaxgenQ']))/bus.Vmin
        #Generous: assumes all power routed through bus, and then some. However it still amounts to an assumption.
        irvar[bus] = model.addVar(obj = 0.0, lb = -ubound, ub = ubound,
                                    name = "ir_%d"%(busf.nodeID))
        bus.irvarind = varcount
        varcount += 1
        ijvar[bus] = model.addVar(obj = 0.0, lb = -ubound, ub = ubound,
                                    name = "ij_%d"%(busf.nodeID))
        bus.ijvarind = varcount
        varcount += 1

    logger.info('Added bus current variables.')
    """

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

    # Save variable data
    # alldata['LP']['cffvar']   = cffvar
    # alldata['LP']['irvar']   = irvar
    # alldata['LP']['ijvar']   = ijvar
    alldata["LP"]["irvar_f"] = irvar_f
    alldata["LP"]["irvar_t"] = irvar_t
    alldata["LP"]["ijvar_f"] = ijvar_f
    alldata["LP"]["ijvar_t"] = ijvar_t
    alldata["LP"]["Pvar_f"] = Pvar_f
    alldata["LP"]["Pvar_t"] = Pvar_t
    alldata["LP"]["Qvar_f"] = Qvar_f
    alldata["LP"]["Qvar_t"] = Qvar_t
    alldata["LP"]["GenPvar"] = GenPvar
    alldata["LP"]["GenQvar"] = GenQvar

    logger.info("Added a total of %d variables." % (varcount))


def lpformulator_iv_create_polar_vars(alldata, model, varcount):
    """
    Creates and adds variables for polar IV formulation to a given Gurobi model
    TODO-Dan This is the same function as in grbformulator_ac.py.
             If this is indeed the same function then please remove 1 of the 2 (or just let me know and I'll do it)

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    model : gurobipy.Model
        Gurobi model to be constructed
    varcount : integer
        Number of variables added to the model so far
    """

    logger = logging.getLogger("OpfLogger")
    vvar = {}
    thetavar = {}
    numbuses = alldata["numbuses"]
    buses = alldata["buses"]
    IDtoCountmap = alldata["IDtoCountmap"]
    numbranches = alldata["numbranches"]
    branches = alldata["branches"]
    newvarcount = 0
    fixtolerance = alldata["fixtolerance"]

    logger.info("  Creating variables for polar formulation.")

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

    logger.info("    Assumption. Phase angle diffs between -pi and pi.")

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

    logger.info("    Added %d new variables to handle polar formulation." % newvarcount)

    varcount += newvarcount


def lpformulator_iv_create_efvars(alldata, model, varcount):
    """
    Creates and adds e, f variables for bilinear IV formulation to a given Gurobi model
    TODO-Dan This is the same function as in grbformulator_ac.py.
             If this is indeed the same function then please remove 1 of the 2 (or just let me know and I'll do it)

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    model : gurobipy.Model
        Gurobi model to be constructed
    varcount : integer
        Number of variables added to the model so far
    """

    logger = logging.getLogger("OpfLogger")
    evar = {}
    fvar = {}
    numbuses = alldata["numbuses"]
    buses = alldata["buses"]
    IDtoCountmap = alldata["IDtoCountmap"]
    efvarcount = 0
    fixtolerance = 1e-05
    if alldata["fixtolerance"] > 0:
        fixtolerance = alldata["fixtolerance"]

    logger.info("  Creating e,f variables.")

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

    logger.info("  Added %d e, f variables." % efvarcount)

    varcount += efvarcount


def lpformulator_iv_create_constraints(alldata, model):
    """
    Creates and adds constraints for IV formulation to a given Gurobi model
    TODO-Dan Is this function very different from the one in grbformulator_ac.py?
             Couldn't we call lpformulator_ac_create_constraints(alldata, model)
             and then add the missing IV constraints?

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    model : gurobipy.Model
        Gurobi model to be constructed
    """

    logger = logging.getLogger("OpfLogger")
    numbuses = alldata["numbuses"]
    buses = alldata["buses"]
    numbranches = alldata["numbranches"]
    branches = alldata["branches"]
    gens = alldata["gens"]
    IDtoCountmap = alldata["IDtoCountmap"]
    evar = alldata["LP"]["evar"]
    fvar = alldata["LP"]["fvar"]
    # cffvar        = alldata['LP']['cffvar']
    # irvar        = alldata['LP']['irvar']
    # ijvar        = alldata['LP']['ijvar']
    irvar_f = alldata["LP"]["irvar_f"]
    irvar_t = alldata["LP"]["irvar_t"]
    ijvar_f = alldata["LP"]["ijvar_f"]
    ijvar_t = alldata["LP"]["ijvar_t"]
    Pvar_f = alldata["LP"]["Pvar_f"]
    Pvar_t = alldata["LP"]["Pvar_t"]
    Qvar_f = alldata["LP"]["Qvar_f"]
    Qvar_t = alldata["LP"]["Qvar_t"]

    GenPvar = alldata["LP"]["GenPvar"]
    GenQvar = alldata["LP"]["GenQvar"]
    lincostvar = alldata["LP"]["lincostvar"]

    logger.info("Creating constraints.")
    logger.info("  Adding cost definition.")

    coeff = [gen.costvector[gen.costdegree - 1] for gen in gens.values()]
    variables = [GenPvar[gen] for gen in gens.values()]
    expr = gp.LinExpr(coeff, variables)
    model.addConstr(expr == lincostvar, name="lincostdef")

    numquadgens = 0
    for gen in gens.values():
        if gen.costdegree >= 2 and gen.costvector[0] > 0 and gen.status:
            numquadgens += 1

    logger.info(
        "    Number of generators with quadratic cost coefficient: %d." % numquadgens
    )

    if numquadgens > 0:
        if alldata["usequadcostvar"]:
            quadcostvar = alldata["LP"]["quadcostvar"]
            logger.info("    Adding quadcost definition constraint.")
            qcost = gp.QuadExpr()
            for gen in gens.values():
                if gen.costdegree == 2 and gen.costvector[0] != 0:
                    qcost.add(gen.costvector[0] * GenPvar[gen] * GenPvar[gen])

            model.addConstr(qcost <= quadcostvar, name="qcostdef")
        else:
            logger.info("    Adding quad cost to objective.")
            model.update()  # Necessary to flush changes in the objective function
            oldobj = model.getObjective()
            newobj = gp.QuadExpr(oldobj)
            for gen in gens.values():
                if gen.costdegree == 2 and gen.costvector[0] != 0:
                    newobj.add(gen.costvector[0] * GenPvar[gen] * GenPvar[gen])

            model.setObjective(newobj, GRB.MINIMIZE)

    if alldata["ivtype"] == "plain":
        logger.info("  Adding branch current definitions.")
        # (g + jb)(e + jf) = g*e - b*f + j( b*e + g*f )
        count = 0
        for j in range(1, 1 + numbranches):
            branch = branches[j]
            f = branch.f
            t = branch.t
            count_of_f = IDtoCountmap[f]
            count_of_t = IDtoCountmap[t]
            busf = buses[count_of_f]
            bust = buses[count_of_t]
            branch.Irfcname = "Irdef_%d_%d_%d" % (j, f, t)
            branch.Irtcname = "Irdef_%d_%d_%d" % (j, t, f)

            branch.Ijfcname = "Ijdef_%d_%d_%d" % (j, f, t)
            branch.Ijtcname = "Ijdef_%d_%d_%d" % (j, t, f)

            if branch.status:
                expr = (
                    branch.Gff * evar[busf]
                    - branch.Bff * fvar[busf]
                    + branch.Gft * evar[bust]
                    - branch.Bft * fvar[bust]
                )
            model.addConstr(expr == irvar_f[branch], name=branch.Irfcname)

            expr = (
                branch.Gtf * evar[busf]
                - branch.Btf * fvar[busf]
                + branch.Gtt * evar[bust]
                - branch.Btt * fvar[bust]
            )
            model.addConstr(expr == irvar_t[branch], name=branch.Irtcname)

            expr = (
                branch.Gff * fvar[busf]
                + branch.Bff * evar[busf]
                + branch.Gft * fvar[bust]
                + branch.Bft * evar[bust]
            )
            model.addConstr(expr == ijvar_f[branch], name=branch.Ijfcname)

            expr = (
                branch.Gtf * fvar[busf]
                + branch.Btf * evar[busf]
                + branch.Gtt * fvar[bust]
                + branch.Btt * evar[bust]
            )
            model.addConstr(expr == ijvar_t[branch], name=branch.Ijtcname)

            count += 4

        logger.info("    %d branch current definitions added." % count)

        logger.info("  Adding plain branch power definitions.")
        # (e + jf)*(Ir - jIj) = e*Ir + f*Ij + j( -e*Ij + f*Ir)
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
            branch.Qfcname = "Qdef_%d_%d_%d" % (j, f, t)
            branch.Qtcname = "Qdef_%d_%d_%d" % (j, t, f)

            if branch.status:
                qexpr = evar[busf] * irvar_f[branch] + fvar[busf] * ijvar_f[branch]
                model.addConstr(qexpr == Pvar_f[branch], name=branch.Pfcname)

                qexpr = evar[bust] * irvar_t[branch] + fvar[bust] * ijvar_t[branch]
                model.addConstr(qexpr == Pvar_t[branch], name=branch.Ptcname)

                qexpr = -evar[busf] * ijvar_f[branch] + fvar[busf] * irvar_f[branch]
                model.addConstr(qexpr == Qvar_f[branch], name=branch.Qfcname)

                qexpr = -evar[bust] * ijvar_t[branch] + fvar[bust] * irvar_t[branch]
                model.addConstr(qexpr == Qvar_t[branch], name=branch.Qtcname)

            count += 4
        logger.info("    %d branch power flow definitions added." % count)
    else:
        logger.info(
            "Aggressive formulation, so will define power flows as functions of voltages."
        )

        logger.info("  Adding aggressive branch power definitions.")
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
            branch.Qfcname = "Qdef_%d_%d_%d" % (j, f, t)
            branch.Qtcname = "Qdef_%d_%d_%d" % (j, t, f)

            if branch.status:
                gkk = branch.Gff
                bkk = branch.Bff
                gkm = branch.Gft
                bkm = branch.Bft

                # (e + jf)*(Ir - jIj) = e*Ir + f*Ij + j( -e*Ij + f*Ir)
                # But Ir = gkk*ek - bkk*fk + gkm*em - bkm*fm,
                # And Ij = bkk*ek + gkk*fk + bkm*em + gkm*fm

                # So real power from k to m = ek*[ gkk*ek - bkk*fk + gkm*em - bkm*fm ]
                #                           + fk*[ bkk*ek + gkk*fk + bkm*em + gkm*fm ] =
                #
                #                   gkk*(ek^2 + fk^2) + gkm*(ek*em + fk*fm) + bkm*(-ek*fm + em*fk)

                qexpr = gkk * (evar[busf] * evar[busf] + fvar[busf] * fvar[busf])
                qexpr += gkm * (evar[busf] * evar[bust] + fvar[busf] * fvar[bust])
                qexpr += bkm * (-evar[busf] * fvar[bust] + fvar[busf] * evar[bust])

                model.addConstr(qexpr == Pvar_f[branch], name=branch.Pfcname)

                # And imag power from k to m = -ek*[ bkk*ek + gkk*fk + bkm*em + gkm*fm ]
                #                            +  fk*[ gkk*ek - bkk*fk + gkm*em - bkm*fm ] =
                #
                #                 -bkk*(ek^2 + fk^2) - bkm*( ek*em + fk*fm) + gkm*(-ek*fm + em*fk)

                qexpr = -bkk * (evar[busf] * evar[busf] + fvar[busf] * fvar[busf])
                qexpr += -bkm * (evar[busf] * evar[bust] + fvar[busf] * fvar[bust])
                qexpr += gkm * (-evar[busf] * fvar[bust] + fvar[busf] * evar[bust])

                model.addConstr(qexpr == Qvar_f[branch], name=branch.Qfcname)

                # now, the reversals

                gmm = branch.Gtt
                bmm = branch.Btt
                gmk = branch.Gtf
                bmk = branch.Btf

                qexpr = gmm * (evar[bust] * evar[bust] + fvar[bust] * fvar[bust])
                qexpr += gmk * (evar[bust] * evar[busf] + fvar[bust] * fvar[busf])
                qexpr += bmk * (-evar[bust] * fvar[busf] + fvar[bust] * evar[busf])

                model.addConstr(qexpr == Pvar_t[branch], name=branch.Ptcname)

                qexpr = -bmm * (evar[bust] * evar[bust] + fvar[bust] * fvar[bust])
                qexpr += -bmk * (evar[bust] * evar[busf] + fvar[bust] * fvar[busf])
                qexpr += gmk * (-evar[bust] * fvar[busf] + fvar[bust] * evar[busf])

                model.addConstr(qexpr == Qvar_t[branch], name=branch.Qtcname)

                count += 4

        logger.info("    %d branch power flow definitions added." % count)

    # Balance constraints

    logger.info("  Adding active power balance constraints.")
    count = 0
    for j in range(1, 1 + numbuses):
        bus = buses[j]
        injpowerexpr = gp.LinExpr()
        for branchid in bus.frombranchids.values():
            injpowerexpr.add(Pvar_f[branches[branchid]])

        for branchid in bus.tobranchids.values():
            injpowerexpr.add(Pvar_t[branches[branchid]])

        genexpr = gp.LinExpr()
        if len(bus.genidsbycount) > 0:
            for genid in bus.genidsbycount:
                gen = gens[genid]
                genexpr.add(GenPvar[gen])

        model.addConstr(
            genexpr == injpowerexpr + bus.Pd, name="PBaldef%d_%d" % (j, bus.nodeID)
        )

        count += 1

    logger.info("    %d active bus power balance constraints added." % count)

    logger.info("  Adding reactive power balance constraints.")
    count = 0
    for j in range(1, 1 + numbuses):
        bus = buses[j]
        injpowerexpr = gp.LinExpr()
        for branchid in bus.frombranchids.values():
            injpowerexpr.add(Qvar_f[branches[branchid]])

        for branchid in bus.tobranchids.values():
            injpowerexpr.add(Qvar_t[branches[branchid]])

        genexpr = gp.LinExpr()
        if len(bus.genidsbycount) > 0:
            for genid in bus.genidsbycount:
                gen = gens[genid]
                genexpr.add(GenQvar[gen])

        model.addConstr(
            genexpr == injpowerexpr + bus.Qd, name="QBaldef%d_%d" % (j, bus.nodeID)
        )

        count += 1

    logger.info("    %d reactive bus power balance constraints added." % count)

    # the next set of constraints is optional
    """
    logger.info("  Adding constraints stating bus current injection = total outgoing current.")
    count = 0
    for j in range(1,1+numbuses):
        bus  = buses[j]
        expr = gp.LinExpr()

        if bus.Gs != 0:
            expr.add(bus.Gs*evar[bus])
        if bus.Bs != 0:
            expr.add(-bus.Bs*fvar[bus])

        if alldata['branchswitching_comp'] == False:
            for branchid in bus.frombranchids.values():
                expr.add(irvar_f[branches[branchid]])

            for branchid in bus.tobranchids.values():
                expr.add(irvar_t[branches[branchid]])
        model.addConstr(expr == irvar[bus], name = "IrBaldef%d_%d"%(j, bus.nodeID))

        count += 1

    for j in range(1,1+numbuses):
        bus  = buses[j]
        expr = gp.LinExpr()

        if bus.Gs != 0:
            expr.add(bus.Gs*fvar[bus])
        if bus.Bs != 0:
            expr.add(bus.Bs*evar[bus])

        if alldata['branchswitching_comp'] == False:
            for branchid in bus.frombranchids.values():
                expr.add(ijvar_f[branches[branchid]])

            for branchid in bus.tobranchids.values():
                expr.add(ijvar_t[branches[branchid]])
        model.addConstr(expr == ijvar[bus], name = "IrBaldef%d_%d"%(j, bus.nodeID))

        count += 1

    logger.info("    %d bus current balance constraints added."%count)
    """

    # Bus voltage limits.

    logger.info("  Adding voltage limits.")

    count = 0
    for j in range(1, 1 + numbuses):
        bus = buses[j]
        model.addConstr(
            evar[bus] ** 2 + fvar[bus] ** 2 <= bus.Vmax * bus.Vmax, name="Vmax_%d" % j
        )
        count += 1
        model.addConstr(
            evar[bus] ** 2 + fvar[bus] ** 2 >= bus.Vmin * bus.Vmin, name="Vmin_%d" % j
        )
        count += 1

    logger.info("    %d bus voltage limit constraints added." % count)

    # Branch limits.
    logger.info("  Adding branch limits.")
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
            # themodel.cbLazy(Pvar_t[branch]*Pvar_t[branch] + Qvar_t[branch]*Qvar_t[branch] <= branch.limit**2)
            model.addConstr(
                Pvar_t[branch] * Pvar_t[branch] + Qvar_t[branch] * Qvar_t[branch]
                <= branch.limit**2,
                name="limit_t_%d_%d_%d" % (j, t, f),
            )
            count += 2

    logger.info("    %d branch limits added." % count)

    # Active loss inequalities.
    if alldata["useactivelossineqs"] == True:
        logger.info("  Adding active loss constraints in weak form.")
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
                    Pvar_f[branch] + Pvar_t[branch] >= 0,
                    name="aLa_%d_%d_%d" % (j, f, t),
                )

                count += 1

        logger.info("    %d active loss inequalities added." % count)

    # break_exit("after linelimits")


def lpformulator_iv_add_polarconstraints(alldata, model):
    """
    Creates and adds constraints for polar represenation of IV formulation to a given Gurobi model
    TODO-Dan This is the same function as in grbformulator_ac.py.
             If this is indeed the same function then please remove 1 of the 2 (or just let me know and I'll do it)

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    model : gurobipy.Model
        Gurobi model to be constructed
    """

    logger = logging.getLogger("OpfLogger")
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

    logger.info("  Adding polar constraints.")

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

    logger.info("    %d polar constraints added." % count)


def lpformulator_iv_add_nonconvexconstraints(alldata, model):
    """
    Creates and adds nonconvex bilinear e, f constraints for AC formulation to a given Gurobi model
    TODO-Dan This is the same function as in grbformulator_ac.py.
             If this is indeed the same function then please remove 1 of the 2 (or just let me know and I'll do it)

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    model : gurobipy.Model
        Gurobi model to be constructed
    """

    logger = logging.getLogger("OpfLogger")
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
    logger.info("  Adding nonconvex e, f, constraints.")

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

    logger.info("    %d nonconvex e, f constraints added." % count)


def lpformulator_iv_examine_solution(alldata, model):
    """
    first version -- crude -- but first, some debugging
    TODO-Dan What is this function and how should it be used?
    Currently it does nothing

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    model : gurobipy.Model
        Constructed Gurobi model
    """

    logger = logging.getLogger("OpfLogger")
    branches = alldata["branches"]
    numbranches = alldata["numbranches"]
    IDtoCountmap = alldata["IDtoCountmap"]

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


                logger.info('j %d f %d Pft %8.4f Qft %8.4f cff %8.4f aL %.2e r %.2e iaL %.3e I2 %.3e aL_r %.3e'%(j,f,Pvar_f[branch].x, Qvar_f[branch].x, cvar[busf].x, aLoss, branch.r, impliedaLoss, actualI2, actualLoss_r))
                #logger.info('          Qft %8.4f Qtf %8.4f rL %.2e x %.2e\n'%(Qvar_f[branch].x, Qvar_t[branch].x, rLoss, branch.x))
                if aLoss < minaloss:
                    minaloss = aLoss
                    minalossind = j



        logger.info('minaloss %.3e at %d'%(minaloss,minalossind))

        """

    logger.info("Done examining solution.\n")


def computebalbounds(alldata, bus):
    """
    Computes active and reactive max and min bus flow balance values

    TODO-Dan There is a function with the same name in grbformulator_iv.
    It looks like the other function does something slightly different.
    Can we combine them to 1 function? Or even delete 1 of the 2

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    bus : Bus
        Input bus
    """

    logger = logging.getLogger("OpfLogger")
    # First let's get max/min generations
    loud = 0  # See below

    baseMVA = alldata["baseMVA"]
    gens = alldata["gens"]
    Pubound = Plbound = 0
    Qubound = Qlbound = 0

    for gencounter in bus.genidsbycount:
        if gens[gencounter].status:
            Pubound += gens[gencounter].Pmax
            Plbound += gens[gencounter].Pmin
            Qubound += gens[gencounter].Qmax
            Qlbound += gens[gencounter].Qmin

        if (
            loud
        ):  # it is worth keeping because a user may want to debug the generator limits that they imposed, which could make a problem infeasible
            # logger.info(" Pubound for %d %f genc %d."%(bus.nodeID, Pubound, gencounter))
            # logger.info(" Plbound for %d %f genc %d."%(bus.nodeID, Plbound, gencounter))
            logger.info(
                " Qubound for %d %f genc %d." % (bus.nodeID, Qubound, gencounter)
            )
            logger.info(
                " Qlbound for %d %f genc %d." % (bus.nodeID, Qlbound, gencounter)
            )

    Pubound -= bus.Pd
    Plbound -= bus.Pd
    Qubound -= bus.Qd
    Qlbound -= bus.Qd

    if bus.nodetype == 4:
        Pubound = Plbound = Qubound = Qlbound = 0

    if loud:  # same as above
        # logger.info(" Pubound for %d final %f."%(bus.nodeID, Pubound))
        # logger.info(" (Pd was %g)\n"%bus.Pd)
        # logger.info(" Plbound for %d final %f."%(bus.nodeID, Plbound))
        logger.info(" Qubound for %d final %f." % (bus.nodeID, Qubound))
        logger.info(" (Qd was %g)" % bus.Qd)
        logger.info(" Qlbound for %d final %f." % (bus.nodeID, Qlbound))

    return Pubound, Plbound, Qubound, Qlbound


def grbderive_xtra_sol_values_fromvoltages(alldata):
    """
    Generates complete solution vectors from input voltages
    TODO-Dan This is the same function as in grbformulator_ac.py.
             If this is indeed the same function then please remove 1 of the 2 (or just let me know and I'll do it)

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    """

    logger = logging.getLogger("OpfLogger")
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

    logger.info("Creating derived solution values from input voltage solution.")

    loud = False

    for j in range(1, numbuses + 1):
        bus = buses[j]
        bus.inpute = bus.inputV * math.cos(bus.inputA_rad)
        bus.inputf = bus.inputV * math.sin(bus.inputA_rad)

    if alldata["use_ef"]:
        evar = alldata["LP"]["evar"]
        fvar = alldata["LP"]["fvar"]
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

    logger.info("Derived e, f, c, s values.")

    Pvar_f = alldata["LP"]["Pvar_f"]
    Pvar_t = alldata["LP"]["Pvar_t"]
    Qvar_f = alldata["LP"]["Qvar_f"]
    Qvar_t = alldata["LP"]["Qvar_t"]

    for j in range(1, 1 + numbranches):
        branch = branches[j]
        row = model.getRow(branch.Pdeffconstr)
        if loud:
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
        logger.info(
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
        logger.info(
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
        logger.info(
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
        logger.info(
            "%s derived at %f." % (Qvar_t[branch].Varname, xbuffer[Qvar_t[branch]])
        )
    # break_exit('branchPQ')

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

    # break_exit('derive')


def lpformulator_iv_strictchecker(alldata, model, spitoutvector):
    """
    Check feasibility of input solution -- report infeasibilities
    TODO-Dan How do we use this function? It says "check input solution against formulation".
             Which input solution? Where does the user provide the input solution? Is it in the voltage file?
    TODO-Dan This is the same function as in grbformulator_ac.py.
             If this is indeed the same function then please remove 1 of the 2 (or just let me know and I'll do it)

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    model : gurobipy.Model
        Gurobi model to be constructed
    spitoutvector: boolean
        # TODO-Dan What does it do?
    """

    logger = logging.getLogger("OpfLogger")
    logger.info("Strict feasibility check for input voltage solution.\n")

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

    logger.info("Direct bus magnitude bounds check. Error issued if infeasibility.")

    for j in range(1, numbuses + 1):
        bus = buses[j]
        if bus.inputV > bus.Vmax:
            logger.error(
                ">>> Error: bus # %d has input voltage %f which is larger than Vmax %f."
                % (j, bus.inputV, bus.Vmax)
            )
            thisviol = bus.inputV - bus.Vmax
            if thisviol > max_violation_value:
                max_violation_string = "bus_" + str(bus.nodeID) + "_Vmax"
                max_violation_value = thisviol
        if bus.inputV < bus.Vmin:
            logger.error(
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

    logger.info("Checked input Vmag values.")

    logger.info("Direct branch limit check. Error issued if infeasible.")

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
            logger.error(
                ">>> Error: branch # %d has 'from' flow magnitude %f which is larger than limit %f."
                % (j, fromvalue, branch.limit)
            )
            logger.error("    branch is ( %d %d )." % (branch.f, branch.t))
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
            logger.error(
                ">>> Error: branch # %d has 'to' flow magnitude %f which is larger than limit %f."
                % (j, tovalue, branch.limit)
            )
            logger.error("    branch is ( %d %d )." % (branch.f, branch.t))
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

        logger.info("e, f values.")

    if alldata["dopolar"]:
        logger.info(
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

        logger.info("Vmag values checked.")

    logger.info("Checking power flow values.")

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

    logger.info("Checking flow balance constraints.")
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

        logger.info(
            "Bus ID %s #%d injection %g mingen %g maxgen %g load %g."
            % (bus.nodeID, j, injection, myPlbound, myPubound, bus.Pd)
        )
        logger.info("   min net generation %g; max %g." % (minnetgen, maxnetgen))

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

        # but also, construct min/max injection at the bus by looking at available generators and load
        myQubound = myQlbound = 0
        for gencounter in bus.genidsbycount:
            if gens[gencounter].status:
                myQubound += gens[gencounter].Qmax
                myQlbound += gens[gencounter].Qmin

        minnetgen = myQlbound - bus.Qd
        maxnetgen = myQubound - bus.Qd

        logger.info(
            "Bus ID %s #%d injection %g mingen %g maxgen %g load %g."
            % (bus.nodeID, j, injection, myQlbound, myQubound, bus.Qd)
        )
        logger.info("   min net generation %g; max %g." % (minnetgen, maxnetgen))

        candmaxviol = alldata["violation"][bus]["Qinjmax"]
        if candmaxviol < alldata["violation"][bus]["Qinjmin"]:
            candmaxviol = -alldata["violation"][bus]["Qinjmin"]
        IQviol[bus] = candmaxviol

    worstboundviol_report(badlbvar, maxlbviol, "LB")
    worstboundviol_report(badubvar, maxubviol, "UB")

    logger.info("\nSummary: Max LB viol %g, Max UB viol %g." % (maxlbviol, maxubviol))
    logger.info(
        "\nSummary: Max violation %g, key: %s."
        % (max_violation_value, max_violation_string)
    )

    # if alldata["dographics"]:
    #    textlist = []
    #    grbgraphical(alldata, "violation", textlist)


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
    loud,  # TODO-Dan Remove loud?
):
    """
    Return bounds infeasibility if setting grbvariable to some value
    maxlbviol, maxubviol, badlb, ubvar are updated if the infeasibility is larger
    TODO-Dan This function definitely looks the same as in grbformulator_ac
             Are these the same function? If yes, please remove 1 of the 2.
             Wouldn't it be enough to use Gurobi attributes BoundVio?

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    model : gurobipy.Model
        Gurobi model to be constructed
    grbvariable : gurobipy.Var
        Gurobi variable of interest
    value : double
        Value of optimization variable of interest
    maxlbviol : double
        Current max lower bound violation over all variables
    maxubviol : double
        Current max upper bound violation over all variables
    badlbvar : gurobipy.Var
        Variable with the biggest lower bound violation
    badubvar : gurobipy.Var
        Variable with the biggest upper bound violation
    max_violation_value : double
        Max overall bound violation # TODO-Dan We can get this via model.BoundVio, see https://www.gurobi.com/documentation/current/refman/attributes.html
    max_violation_string : string
        Name of variable with largest overall bound violation # TODO-Dan We can get this via model.BoundVioIndex, see https://www.gurobi.com/documentation/current/refman/attributes.html
    """

    logger = logging.getLogger("OpfLogger")
    ub = grbvariable.ub
    lb = grbvariable.lb

    if loud:
        logger.info(
            "%s =  %.16e  [ LB %.4e  UB %.4e ]" % (grbvariable.varname, value, lb, ub)
        )

    lbviol = lb - value
    ubviol = value - ub

    if lbviol > 0:
        logger.info(
            "LBVIOL %s  LB %.6e  x %.16e  UB %.6e"
            % (grbvariable.varname, lb, value, ub)
        )

    if lbviol > maxlbviol:
        logger.info(" -> incumbent MAX LBVIOL")
        maxlbviol = lb - value
        badlbvar = grbvariable

    if ubviol > 0:
        logger.info(
            "UBVIOL %s  LB %.6e  x %.16e  UB %.6e"
            % (grbvariable.varname, lb, value, ub)
        )

    if ubviol > maxubviol:
        logger.info(" -> incumbent MAX UBVIOL")
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
    """
    Report the variable with largest bound violation
    TODO-Dan This function definitely looks the same as in grbformulator_ac
             Are these the same function? If yes, please remove 1 of the 2.

    Parameters
    ----------
    badvar : gurobipy.Var
        Gurobi variable with largest bound violation
    maxviol : double
        Value of violation
    boundtype : string
        States whether it's a lower or an upper bound
    """

    logger = logging.getLogger("OpfLogger")
    if badvar != None:
        logger.info(
            "Worst %s bound violation by %s viol %g."
            % (boundtype, badvar.Varname, maxviol)
        )
    else:
        logger.info("No %s bound violations." % boundtype)
