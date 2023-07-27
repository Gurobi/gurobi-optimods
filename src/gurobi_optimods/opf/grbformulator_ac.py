import logging
import math

import gurobipy as gp
from gurobipy import GRB

from gurobi_optimods.opf.grbformulator_common import set_gencost_objective

logger = logging.getLogger(__name__)


def lpformulator_ac_body(alldata, model):
    """
    Adds variables and constraints for AC formulation to a given Gurobi model

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param model: Gurobi model to be constructed
    :type model: :class: `gurobipy.Model`
    """

    lpformulator_ac_create_vars(alldata, model)
    set_gencost_objective(alldata, model)
    lpformulator_ac_create_constraints(alldata, model)


def lpformulator_ac_create_vars(alldata, model):
    """
    Creates and adds variables for AC formulation to a given Gurobi model

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param model: Gurobi model to be constructed
    :type model: :class: `gurobipy.Model`
    """

    logger.info("Creating variables.")

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

    for j in range(1, numbuses + 1):
        bus = buses[j]
        # First, injection variables
        maxprod = bus.Vmax * bus.Vmax
        minprod = bus.Vmin * bus.Vmin
        ubound = maxprod
        lbound = minprod

        # Check whether to fix variable and we have a input voltage solution
        if alldata["fixcs"] and bus.inputvoltage:
            lbound = bus.inputV * bus.inputV - fixtolerance
            ubound = bus.inputV * bus.inputV + fixtolerance

        cvar[bus] = model.addVar(
            obj=0.0, lb=lbound, ub=ubound, name="c_%d_%d" % (bus.nodeID, bus.nodeID)
        )

        # csdefslacks to be done
        Plbound = Qlbound = -GRB.INFINITY
        Pubound = Qubound = GRB.INFINITY
        Pubound, Plbound, Qubound, Qlbound = computebalbounds(alldata, bus)

        Pinjvar[bus] = model.addVar(
            obj=0.0, lb=Plbound, ub=Pubound, name="IP_%d" % bus.nodeID
        )
        # comment: Pinjvar is the variable modeling total active power injected by bus j into the branches incident with j

        Qinjvar[bus] = model.addVar(
            obj=0.0, lb=Qlbound, ub=Qubound, name="IQ_%d" % bus.nodeID
        )
        # comment: Qinjvar is the variable modeling total reactive power injected by bus j into the branches incident with j

        # Next, generator variables in a bus
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

            lower = gen.Qmin * gen.status
            upper = gen.Qmax * gen.status
            if bus.nodetype == 3:
                lower = -GRB.INFINITY
                upper = GRB.INFINITY

            GenQvar[gen] = model.addVar(
                obj=0.0, lb=lower, ub=upper, name="GQ_%d_%d" % (gen.count, gen.nodeID)
            )

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
            logger.error(
                f" --- Broken assumption 1: branch j {j} f {f} t {t} minanglerad {branch.minangle_rad} maxanglerad {branch.maxangle_rad}."
            )
            raise ValueError(
                f"Broken assumption 1: branch j {j} f {f} t {t} minanglerad {branch.minangle_rad} maxanglerad {branch.maxangle_rad}."
            )

        ubound = maxprod
        lbound = -maxprod
        maxanglerad = branch.maxangle_rad
        minanglerad = branch.minangle_rad

        # Cosine
        if maxanglerad <= 0.5 * math.pi:
            # In this case minangle <= 0
            if minanglerad >= -0.5 * math.pi:
                lbound = minprod * min(math.cos(maxanglerad), math.cos(minanglerad))
            elif minanglerad >= -math.pi:
                lbound = maxprod * math.cos(minanglerad)  # Which is negative
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

        Pvar_t[branch] = model.addVar(
            obj=0.0,
            lb=lbound,
            ub=ubound,
            name="P_%d_%d_%d" % (j, bust.nodeID, busf.nodeID),
        )

        Qvar_f[branch] = model.addVar(
            obj=0.0,
            lb=lbound,
            ub=ubound,
            name="Q_%d_%d_%d" % (j, busf.nodeID, bust.nodeID),
        )

        Qvar_t[branch] = model.addVar(
            obj=0.0,
            lb=lbound,
            ub=ubound,
            name="Q_%d_%d_%d" % (j, bust.nodeID, busf.nodeID),
        )

        if alldata["branchswitching_mip"]:
            twinPvar_f[branch] = model.addVar(
                obj=0.0,
                lb=lbound,
                ub=ubound,
                name="twinP_%d_%d_%d" % (j, busf.nodeID, bust.nodeID),
            )

            twinPvar_t[branch] = model.addVar(
                obj=0.0,
                lb=lbound,
                ub=ubound,
                name="twinP_%d_%d_%d" % (j, bust.nodeID, busf.nodeID),
            )

            twinQvar_f[branch] = model.addVar(
                obj=0.0,
                lb=lbound,
                ub=ubound,
                name="twinQ_%d_%d_%d" % (j, busf.nodeID, bust.nodeID),
            )

            twinQvar_t[branch] = model.addVar(
                obj=0.0,
                lb=lbound,
                ub=ubound,
                name="twinQ_%d_%d_%d" % (j, bust.nodeID, busf.nodeID),
            )

    zvar = {}
    if alldata["branchswitching_mip"] or alldata["branchswitching_comp"]:
        logger.info("Adding branch switching variables.")
        for j in range(1, 1 + numbranches):
            branch = branches[j]
            f = branch.f
            t = branch.t
            zvar[branch] = model.addVar(
                obj=0.0,
                vtype=GRB.BINARY,
                name="z_%d_%d_%d" % (j, f, t),
            )

    # Powerflow variables
    if alldata["use_ef"]:
        lpformulator_ac_create_efvars(alldata, model)

    if alldata["dopolar"]:
        lpformulator_ac_create_polar_vars(alldata, model)

    # Save variable data
    alldata["LP"]["cvar"] = cvar
    alldata["LP"]["svar"] = svar
    alldata["LP"][
        "Pvar_f"
    ] = Pvar_f  # AC branch real power injected into "from" end of branch
    alldata["LP"][
        "Pvar_t"
    ] = Pvar_t  # AC branch real power injected into "to" end of branch
    alldata["LP"][
        "Qvar_f"
    ] = Qvar_f  # AC branch reactive power injected into "from" end of branch
    alldata["LP"][
        "Qvar_t"
    ] = Qvar_t  # AC branch reactive power injected into "to" end of branch
    alldata["MIP"]["zvar"] = zvar
    if alldata["branchswitching_mip"]:
        alldata["LP"]["twinPvar_f"] = twinPvar_f
        alldata["LP"]["twinPvar_t"] = twinPvar_t
        alldata["LP"]["twinQvar_f"] = twinQvar_f
        alldata["LP"]["twinQvar_t"] = twinQvar_t

    alldata["LP"]["GenPvar"] = GenPvar  # AC generator real power injections
    alldata["LP"]["GenQvar"] = GenQvar  # AC generator reactive power injections
    alldata["LP"][
        "Pinjvar"
    ] = Pinjvar  # Pinjvar is the variable modeling total active power injected by bus j into the branches incident with j
    alldata["LP"][
        "Qinjvar"
    ] = Qinjvar  # Qinjvar is the variable modeling total reactive power injected by bus j into the branches incident with j


def lpformulator_ac_create_polar_vars(alldata, model):
    """
    Creates and adds variables for polar AC formulation to a given Gurobi model

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param: Gurobi model to be constructed
    :type model: :class: `gurobipy.Model`
    """

    numbuses = alldata["numbuses"]
    buses = alldata["buses"]
    IDtoCountmap = alldata["IDtoCountmap"]
    numbranches = alldata["numbranches"]
    branches = alldata["branches"]
    fixtolerance = alldata["fixtolerance"]
    newvarcount = 0

    logger.info("  Creating variables for polar formulation.")

    vvar = {}
    thetavar = {}
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
    alldata["LP"]["vvar"] = vvar  # voltage magnitude of bus for polar case
    alldata["LP"]["thetavar"] = thetavar  # voltage angle of bus for polar case
    alldata["LP"]["cosvar"] = cosvar
    alldata["LP"]["sinvar"] = sinvar
    alldata["LP"]["thetaftvar"] = thetaftvar
    alldata["LP"]["vfvtvar"] = vfvtvar

    logger.info(f"    Added {newvarcount} new variables to handle polar formulation.")


def lpformulator_ac_create_efvars(alldata, model):
    """
    Creates and adds e, f variables for bilinear AC formulation to a given Gurobi model

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param model: Gurobi model to be constructed
    :type model: :class: `gurobipy.Model`
    """

    numbuses = alldata["numbuses"]
    buses = alldata["buses"]
    alldata["IDtoCountmap"]
    efvarcount = 0

    logger.info("  Creating e, f variables.")

    evar = {}
    fvar = {}
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

        if alldata["usemaxdispersion"]:
            lbound = 0
            ubound = bus.Vmax * math.sin(alldata["maxdispersion_rad"])
        elif j == alldata["refbus"]:
            ubound = lbound = 0

        fvar[bus] = model.addVar(
            obj=0.0, lb=lbound, ub=ubound, name="f_%d" % bus.nodeID
        )
        efvarcount += 2

    # Save e, f variables data
    alldata["LP"]["evar"] = evar
    alldata["LP"]["fvar"] = fvar  # e^2 + f^2 = voltage magnitude

    logger.info(f"  Added {efvarcount} e, f variables.")


def lpformulator_ac_create_constraints(alldata, model):
    """
    Creates and adds constraints for AC formulation to a given Gurobi model

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param: Gurobi model to be constructed
    :type model: :class: `gurobipy.Model`
    """

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
    zvar = alldata["MIP"]["zvar"]

    logger.info("Creating constraints.")
    logger.info("  Adding cost definition.")

    # Active PF defs
    logger.info("  Adding active power flow definitions.")
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

        if not branch.status:  # out of operation
            branch.Pdeffconstr = model.addConstr(
                Pvar_f[branch] == 0, name=branch.Pfcname
            )
            branch.Pdeftconstr = model.addConstr(
                Pvar_t[branch] == 0, name=branch.Ptcname
            )
            continue

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

    logger.info(f"    {count} active power flow definitions added.")

    # Reactive PF defs
    logger.info("  Adding reactive power flow definitions.")
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

        if not branch.status:  # out of operation
            branch.Qdeffconstr = model.addConstr(
                Qvar_f[branch] == 0, name=branch.Qfcname
            )
            branch.Qdeftconstr = model.addConstr(
                Qvar_t[branch] == 0, name=branch.Qtcname
            )
            continue

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

    logger.info(f"    {count} reactive power flow definitions added.")

    # Balance constraints
    logger.info(
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

        if alldata["branchswitching_comp"] is False:
            for branchid in bus.frombranchids.values():
                expr.add(Pvar_f[branches[branchid]])

            for branchid in bus.tobranchids.values():
                expr.add(Pvar_t[branches[branchid]])
        else:
            doquad = True
            # product of binary and powerflow
            for branchid in bus.frombranchids.values():
                qexpr.add(zvar[branches[branchid]] * Pvar_f[branches[branchid]])

            for branchid in bus.tobranchids.values():
                qexpr.add(zvar[branches[branchid]] * Pvar_t[branches[branchid]])

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

        if alldata["branchswitching_comp"] is False:
            for branchid in bus.frombranchids.values():
                expr.add(Qvar_f[branches[branchid]])

            for branchid in bus.tobranchids.values():
                expr.add(Qvar_t[branches[branchid]])

            model.addConstr(expr == Qinjvar[bus], name="QBaldef%d_%d" % (j, bus.nodeID))
        else:
            # product of binary and powerflow
            qexpr = gp.QuadExpr()

            for branchid in bus.frombranchids.values():
                qexpr.add(zvar[branches[branchid]] * Qvar_f[branches[branchid]])

            for branchid in bus.tobranchids.values():
                qexpr.add(zvar[branches[branchid]] * Qvar_t[branches[branchid]])

            model.addConstr(
                expr + qexpr == Qinjvar[bus], name="QBaldef%d_%d" % (j, bus.nodeID)
            )
        count += 1

    logger.info(f"    {count} balance constraints added.")

    # Injection defs
    logger.info("  Adding injection definition constraints.")
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

    logger.info(f"    {count} injection definition constraints added.")

    # Branch limits
    logger.info("  Adding branch limits.")
    count = 0
    for j in range(1, 1 + numbranches):
        branch = branches[j]

        if not branch.status or branch.unboundedlimit is True:
            continue

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

    logger.info(f"    {count} branch limits added.")

    # JABR.
    if alldata["skipjabr"] is False:
        logger.info("  Adding Jabr constraints.")
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

        logger.info(f"    {count} Jabr constraints added.")
    else:
        logger.info("  Skipping Jabr inequalities.")

    # Active loss constraints.
    if alldata["useactivelossineqs"] is True:
        logger.info("  Adding active loss constraints.")
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

                    # Reactive, also
                    if -branch.Bff > 0 or -branch.Btt > 0:
                        alpha = min(-branch.Bff, -branch.Btt)
                        beta = branch.Bft + branch.Btf
                        delta = 2 * alpha - beta
                        gamma = -branch.Gft + branch.Gtf

                        if delta >= 0:
                            model.addConstr(
                                Qvar_f[branch] + Qvar_t[branch] >= gamma * svar[branch],
                                name="rLa_%d_%d_%d" % (j, f, t),
                            )

                            count += 1

        logger.info(f"    {count} active loss constraints added.")
    else:
        logger.info("  Skipping active loss inequalities.")

    # Polar representation
    if alldata["dopolar"]:
        lpformulator_ac_add_polarconstraints(alldata, model)

    # nonconvex e, f representation
    if alldata["use_ef"]:
        lpformulator_ac_add_nonconvexconstraints(alldata, model)

    if alldata["branchswitching_mip"] or alldata["branchswitching_comp"]:
        expr = gp.LinExpr()
        N = math.floor(
            numbranches * alldata["minactivebranches"]
        )  # <<<<<<---- here is the heuristic lower bound
        logger.info(f"In bound_zs constraint, N = {N}.")
        for j in range(1, 1 + numbranches):
            branch = branches[j]
            expr.add(zvar[branch])
        model.addConstr(expr >= N, name="sumz_lower_heuristic_bound")


def lpformulator_ac_add_polarconstraints(alldata, model):
    """
    Creates and adds constraints for polar represenation of AC formulation to a given Gurobi model

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param model: Gurobi model to be constructed
    :type model: :class: `gurobipy.Model`
    """

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

        model.addGenConstrCos(
            thetaftvar[branch], cosvar[branch], name="cosdef_%d_%d_%d" % (j, f, t)
        )
        model.addGenConstrSin(
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

    logger.info(f"    {count} polar constraints added.")


def lpformulator_ac_add_nonconvexconstraints(alldata, model):
    """
    Creates and adds nonconvex bilinear e, f constraints for AC formulation to a given Gurobi model

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param model: Gurobi model to be constructed
    :type model: :class: `gurobipy.Model`
    """

    buses = alldata["buses"]
    numbuses = alldata["numbuses"]
    branches = alldata["branches"]
    numbranches = alldata["numbranches"]
    evar = alldata["LP"]["evar"]
    fvar = alldata["LP"]["fvar"]  #
    # cvar[km] = (voltage mag at k) * (voltage mag at m) * cos(thetak - thetam)
    # Compute voltage magnitudes for k and m
    # Take arccos to get thetak - thetam = some value
    # Set arbitrarily theta1 = 0 and solve all others
    cvar = alldata["LP"]["cvar"]  # Square of the voltage magnitude e^2 + f^2
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

    logger.info(f"    {count} nonconvex e, f constraints added.")


def computebalbounds(alldata, bus):
    """
    Computes active and reactive max and min bus flow balance values

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param bus: Input bus
    :type bus: :class: `Bus`

    :return: Values for active and reactive bounds
    :rtype: float, float, float, float
    """

    # First let's get max/min generations
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
                upper = GRB.INFINITY
                lower = -GRB.INFINITY

            Qubound += upper
            Qlbound += lower

    Pubound -= bus.Pd
    Plbound -= bus.Pd
    Qubound -= bus.Qd
    Qlbound -= bus.Qd

    if bus.nodetype == 4:
        Pubound = Plbound = Qubound = Qlbound = 0

    return Pubound, Plbound, Qubound, Qlbound


def grbderive_xtra_sol_values_from_voltages(alldata, model):
    """
    Computes complete solution vectors from input voltages
    required for violation computation

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param model: Gurobi model
    :type model: :class: `gurobipy.Model`
    """

    model.update()

    xbuffer = alldata["LP"]["xbuffer"] = {}  # dictionary to store solution values
    numbuses = alldata["numbuses"]
    buses = alldata["buses"]
    IDtoCountmap = alldata["IDtoCountmap"]
    branches = alldata["branches"]
    numbranches = alldata["numbranches"]
    cvar = alldata["LP"]["cvar"]
    svar = alldata["LP"]["svar"]

    for j in range(1, numbuses + 1):
        bus = buses[j]
        bus.inpute = bus.inputV * math.cos(bus.inputA_rad)
        bus.inputf = bus.inputV * math.sin(bus.inputA_rad)

    # Compute xbuffer
    if alldata["use_ef"]:
        evar = alldata["LP"]["evar"]
        fvar = alldata["LP"]["fvar"]
        for j in range(1, numbuses + 1):
            bus = buses[j]
            xbuffer[evar[bus]] = bus.inpute
            xbuffer[fvar[bus]] = bus.inputf
        logger.info("Derived e, f values.")

    if alldata["dopolar"] is False:
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
        # Note: we may not need all of these
        vfvtvar = alldata["LP"]["vfvtvar"]
        thetaftvar = alldata["LP"]["thetaftvar"]
        cosvar = alldata["LP"]["cosvar"]
        sinvar = alldata["LP"]["sinvar"]
        vvar = alldata["LP"]["vvar"]
        thetavar = alldata["LP"]["thetavar"]

        for j in range(1, numbuses + 1):
            bus = buses[j]
            xbuffer[cvar[bus]] = bus.inputV * bus.inputV
            xbuffer[vvar[bus]] = bus.inputV
            xbuffer[thetavar[bus]] = bus.inputA_rad

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

    logger.info("Derived c, s values.")

    Pvar_f = alldata["LP"]["Pvar_f"]
    Pvar_t = alldata["LP"]["Pvar_t"]
    Qvar_f = alldata["LP"]["Qvar_f"]
    Qvar_t = alldata["LP"]["Qvar_t"]

    for j in range(1, 1 + numbranches):
        branch = branches[j]

        for item in [
            (branch.Pdeffconstr, Pvar_f[branch]),
            (branch.Pdeftconstr, Pvar_t[branch]),
            (branch.Qdeffconstr, Qvar_f[branch]),
            (branch.Qdeftconstr, Qvar_t[branch]),
        ]:
            constr = item[0]
            var = item[1]

            row = model.getRow(constr)
            sum = -constr.RHS
            leadcoeff = 0
            for i in range(row.size()):
                v = row.getVar(i)
                coeff = row.getCoeff(i)
                if v.Varname != var.Varname:
                    sum += coeff * xbuffer[v]
                else:
                    leadcoeff = coeff

            xbuffer[var] = -sum / leadcoeff
            # leadcoeff should be +1 or -1

    logger.info("Derived (re)active power flows.")

    # Next, power flow injections
    Pinjvar = alldata["LP"]["Pinjvar"]
    Qinjvar = alldata["LP"]["Qinjvar"]
    for j in range(1, 1 + numbuses):
        bus = buses[j]

        injectionP = 0
        injectionQ = 0
        for branchid in bus.frombranchids.values():
            branch = branches[branchid]
            injectionP += xbuffer[Pvar_f[branch]]
            injectionQ += xbuffer[Qvar_f[branch]]

        for branchid in bus.tobranchids.values():
            branch = branches[branchid]
            injectionP += xbuffer[Pvar_t[branch]]
            injectionQ += xbuffer[Qvar_t[branch]]

        xbuffer[Pinjvar[bus]] = injectionP
        xbuffer[Qinjvar[bus]] = injectionQ

    logger.info("Derived (re)active power flow injections.")


def lpformulator_ac_strictchecker(alldata, model):
    """
    Check feasibility of input solution -- report infeasibilities

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param model: Gurobi model
    :type model: :class: `gurobipy.Model`
    """

    # Derive additional violation information before proceeding
    grbderive_xtra_sol_values_from_voltages(alldata, model)

    buses = alldata["buses"]
    numbuses = alldata["numbuses"]
    branches = alldata["branches"]
    numbranches = alldata["numbranches"]
    gens = alldata["gens"]

    max_violation_string = None
    max_violation_value = 0

    alldata["violation"] = {}  # Dictionary to keep track of violations
    Vmagviol = alldata["violation"]["Vmagviol"] = {}
    IPviol = alldata["violation"]["IPviol"] = {}
    IQviol = alldata["violation"]["IQviol"] = {}
    alldata["violation"]["branchlimit"] = {}
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

    logger.info("Direct bus magnitude bounds check. Warning issued if violated.")

    for j in range(1, numbuses + 1):
        bus = buses[j]
        if bus.inputV > bus.Vmax:
            logger.warning(
                f">>> Warning: bus # {j} has input voltage {bus.inputV} which is larger than Vmax {bus.Vmax}."
            )
            thisviol = bus.inputV - bus.Vmax
            if thisviol > max_violation_value:
                max_violation_string = "bus_" + str(bus.nodeID) + "_Vmax"
                max_violation_value = thisviol
        if bus.inputV < bus.Vmin:
            logger.warning(
                f">>> Warning: bus # {j} has input voltage {bus.inputV} which is smaller than Vmin {bus.Vmin}."
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

    logger.info("Checked input Vm values.")

    logger.info("Direct branch limit check. Warning issued if violated.")

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
            logger.warning(
                f">>> Warning: branch # {j} has 'from' flow magnitude {fromvalue} which is larger than limit {branch.limit}."
            )
            logger.warning(f"    branch is ( {branch.f} {branch.t} ).")
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
            logger.warning(
                f">>> Warning: branch # {j} has 'to' flow magnitude {tovalue} which is larger than limit {branch.limit}."
            )
            logger.warning(f"    branch is ( {branch.f} {branch.t} ).")
            thisviol = toviol
            if thisviol > max_violation_value:
                max_violation_string = "branch_" + str(j) + "_to"
                max_violation_value = thisviol
        alldata["violation"]["branchlimit"][branch] = max(fromviol, toviol)

    logger.info("Checked branch limits.")

    if alldata["use_ef"]:
        logger.info("Checking e, f values.")
        for j in range(1, numbuses + 1):
            bus = buses[j]
            for var in [(evar[bus], bus.inpute), (fvar[bus], bus.inputf)]:
                output_string = f"auxiliary variable {var[0].VarName} used for rectangular formulation defined for bus ID {bus.nodeID}"
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
                    var[0],
                    var[1],
                    maxlbviol,
                    maxubviol,
                    badlbvar,
                    badubvar,
                    max_violation_value,
                    max_violation_string,
                    output_string,
                )

            alldata["LP"]["xbuffer"][evar[bus]] = bus.inpute
            alldata["LP"]["xbuffer"][fvar[bus]] = bus.inputf

        logger.info("e, f values checked.")

    if alldata["dopolar"]:
        logger.info(
            "Checking polar quantities. Note: bounds shifted by input solution."
        )  # which will repeat the Vmag value check, so ...
        vvar = alldata["LP"]["vvar"]
        for j in range(1, numbuses + 1):
            bus = buses[j]
            output_string = f"variable {vvar[bus].VarName} defining voltage magnitude of bus ID {bus.nodeID}"
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
                output_string,
            )

        logger.info("Polar quantities checked.")

    logger.info("Checking power flow values.")

    for j in range(1, 1 + numbranches):
        branch = branches[j]
        for var in [Pvar_f[branch], Pvar_t[branch], Qvar_f[branch], Qvar_t[branch]]:
            power = (
                "real power injected"
                if "P" in var.VarName
                else "reactive power injected"
            )
            direction = "from" if "_f" in var.VarName else "to"
            output_string = f"variable {var.VarName} defining {power} into '{direction}' end of branch ({branch.f},{branch.t})"
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
                var,
                xbuffer[var],
                maxlbviol,
                maxubviol,
                badlbvar,
                badubvar,
                max_violation_value,
                max_violation_string,
                output_string,
            )

    logger.info("Power flow values checked.")

    logger.info("Checking flow balance constraints.")
    Pinjvar = alldata["LP"]["Pinjvar"]
    Qinjvar = alldata["LP"]["Qinjvar"]

    # P variables first
    for j in range(1, 1 + numbuses):
        bus = buses[j]
        varinj = Pinjvar[bus]
        output_string = f"variable {varinj.VarName} defining total real power injected by bus ID {bus.nodeID} into incident branches"
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
            varinj,
            xbuffer[varinj],
            maxlbviol,
            maxubviol,
            badlbvar,
            badubvar,
            max_violation_value,
            max_violation_string,
            output_string,
        )

        alldata["violation"][bus]["Pinjmax"] = max(xbuffer[varinj] - varinj.ub, 0)
        alldata["violation"][bus]["Pinjmin"] = max(varinj.lb - xbuffer[varinj], 0)

        injectionP = 0
        for branchid in bus.frombranchids.values():
            branch = branches[branchid]
            varf = Pvar_f[branch]
            injectionP += xbuffer[varf]

        for branchid in bus.tobranchids.values():
            branch = branches[branchid]
            varf = Pvar_f[branch]
            injectionP += xbuffer[varf]

        # at this point, injectionP is the sum of P flows out of the bus

        # Construct min/max injection at the bus by looking at available generators and load
        myPubound = myPlbound = 0
        for gencounter in bus.genidsbycount:
            if gens[gencounter].status:
                myPubound += gens[gencounter].Pmax
                myPlbound += gens[gencounter].Pmin

        myPlbound - bus.Pd
        myPubound - bus.Pd
        # minnetgen, maxnetgen are (resp.) min and max net generation at the given bus.

        # compute candidate

        candmaxviol = alldata["violation"][bus]["Pinjmax"]
        if candmaxviol < alldata["violation"][bus]["Pinjmin"]:
            candmaxviol = -alldata["violation"][bus]["Pinjmin"]
        IPviol[bus] = candmaxviol

    # Q variables second
    for j in range(1, 1 + numbuses):
        bus = buses[j]
        varinj = Qinjvar[bus]
        varf = Qvar_f[branch]
        output_string = f"variable {varinj.VarName} defining total reactive power injected by bus ID {bus.nodeID} into incident branches"
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
            varinj,
            xbuffer[varinj],
            maxlbviol,
            maxubviol,
            badlbvar,
            badubvar,
            max_violation_value,
            max_violation_string,
            output_string,
        )

        alldata["violation"][bus]["Qinjmax"] = max(xbuffer[varinj] - varinj.ub, 0)
        alldata["violation"][bus]["Qinjmin"] = max(varinj.lb - xbuffer[varinj], 0)

        injectionQ = 0
        for branchid in bus.frombranchids.values():
            branch = branches[branchid]
            injectionQ += xbuffer[varf]

        for branchid in bus.tobranchids.values():
            branch = branches[branchid]
            injectionQ += xbuffer[varf]

        # but also, construct min/max injection at the bus by looking at available generators and load
        myQubound = myQlbound = 0
        for gencounter in bus.genidsbycount:
            if gens[gencounter].status:
                myQubound += gens[gencounter].Qmax
                myQlbound += gens[gencounter].Qmin

        myQlbound - bus.Qd
        myQubound - bus.Qd

        candmaxviol = alldata["violation"][bus]["Qinjmax"]
        if candmaxviol < alldata["violation"][bus]["Qinjmin"]:
            candmaxviol = -alldata["violation"][bus]["Qinjmin"]
        IQviol[bus] = candmaxviol

    logger.info("Flow balance constraints checked.")

    worstboundviol_report(badlbvar, maxlbviol, "LB")
    worstboundviol_report(badubvar, maxubviol, "UB")

    logger.info(f"\nSummary: Max LB viol {maxlbviol:.4e}, Max UB viol {maxubviol:.4e}.")
    logger.info(
        f"         Max overall violation {max_violation_value:.4e}, key: {max_violation_string}."
    )


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
    output_string="",
):
    """
    Returns bounds infeasibility if setting grbvariable to some value.
    Inputs maxlbviol, maxubviol, badlb, ubvar are updated if the infeasibility is larger

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param model: Gurobi model
    :type model: :class: `gurobipy.Model`
    :param grbvariable: Gurobi variable of interest
    :type grbvariable: :class: `gurobipy.Var`
    :param value: Value of optimization variable of interest
    :type value: float
    :param maxlbviol: Current max lower bound violation over all variables
    :type maxlbviol: float
    :param maxubviol: Current max upper bound violation over all variables
    :type maxubviol: float
    :param badlbvar: Variable with the biggest lower bound violation
    :type badlbvar: :class: `gurobipy.Var`
    :param badubvar: Variable with the biggest upper bound violation
    :type badubvar: :class: `gurobipy.Var`
    :param max_violation_value: Max overall bound violation
    :type max_violation_value: float
    :param max_violation_string: Name of variable with largest overall bound violation
    :type max_violation_string: str
    """

    ub = grbvariable.ub
    lb = grbvariable.lb

    if False:  # For debugging only
        logger.info(f"{grbvariable.varname} =  {value}  [ LB {lb}  UB {ub} ]")

    lbviol = lb - value
    ubviol = value - ub

    # lbviol, ubviol: violation of lower and upper bound for given variable

    if lbviol > 0:
        logger.info(
            f"LB violation for {output_string} LB {lb:<14.4e}  x {value:<14.4e}  UB {ub:<14.4e}"
        )

    if lbviol > maxlbviol:
        maxlbviol = lb - value
        badlbvar = grbvariable

    if ubviol > 0:
        logger.info(
            f"UB violation for {output_string} LB {lb:<14.4e}  x {value:<14.4e}  UB {ub:<14.4e}"
        )

    if ubviol > maxubviol:
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
    Reports the variable with largest bound violation

    :param badvar: Gurobi variable with largest bound violation
    :type badvar: :class: `gurobipy.Var`
    :param maxviol: Value of violation
    :type maxviol: float
    :param boundtype: States whether it's a lower or an upper bound
    :type boundtype: str
    """

    if badvar is not None:
        logger.info(
            f"Worst {boundtype} bound violation by {badvar.Varname} viol {maxviol}."
        )
    else:
        logger.info(f"No {boundtype} bound violations.")
