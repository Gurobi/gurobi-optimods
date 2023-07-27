import logging
import math

import gurobipy as gp
from gurobipy import GRB

from gurobi_optimods.opf.grbformulator_common import set_gencost_objective

logger = logging.getLogger(__name__)


def lpformulator_ac_body(alldata, model):
    """Add variables and constraints for AC formulation to the given model"""

    lpformulator_ac_create_vars(alldata, model)
    set_gencost_objective(alldata, model)
    lpformulator_ac_create_constraints(alldata, model)


def lpformulator_ac_create_vars(alldata, model):
    """Add variables for AC formulation to a given Gurobi model"""

    fixtolerance = alldata["fixtolerance"]
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

    for j, bus in buses.items():
        # Injection variables
        if alldata["fixcs"] and bus.inputvoltage:
            # We have a input voltage solution: fix variables within tolerance
            lbound = bus.inputV * bus.inputV - fixtolerance
            ubound = bus.inputV * bus.inputV + fixtolerance
        else:
            ubound = bus.Vmax * bus.Vmax
            lbound = bus.Vmin * bus.Vmin
        cvar[bus] = model.addVar(
            lb=lbound, ub=ubound, name=f"c_{bus.nodeID}_{bus.nodeID}"
        )

        # TODO csdefslacks?
        Pubound, Plbound, Qubound, Qlbound = computebalbounds(alldata, bus)

        # Pinjvar is the variable modeling total active power injected by bus j
        # into the branches incident with j.
        Pinjvar[bus] = model.addVar(lb=Plbound, ub=Pubound, name=f"IP_{bus.nodeID}")

        # Qinjvar is the variable modeling total reactive power injected by bus
        # j into the branches incident with j.
        Qinjvar[bus] = model.addVar(lb=Qlbound, ub=Qubound, name=f"IQ_{bus.nodeID}")

        # Generator variables in a bus (fixed to zero if gen.status == 0)
        for genid in bus.genidsbycount:
            gen = gens[genid]

            # ignoring slack bus
            # if bus.nodetype == 3:
            #     upper = GRB.INFINITY
            #     lower = -GRB.INFINITY
            lower = gen.Pmin * gen.status
            upper = gen.Pmax * gen.status
            GenPvar[gen] = model.addVar(
                lb=lower, ub=upper, name=f"GP_{gen.count}_{gen.nodeID}"
            )

            lower = gen.Qmin * gen.status
            upper = gen.Qmax * gen.status
            if bus.nodetype == 3:
                lower = -GRB.INFINITY
                upper = GRB.INFINITY
            GenQvar[gen] = model.addVar(
                lb=lower, ub=upper, name=f"GQ_{gen.count}_{gen.nodeID}"
            )

    # Branch related variables
    branches = alldata["branches"]
    for j, branch in branches.items():
        busf = buses[IDtoCountmap[branch.f]]
        bust = buses[IDtoCountmap[branch.t]]
        maxprod = busf.Vmax * bust.Vmax
        minprod = busf.Vmin * bust.Vmin

        # Assumption 1.  zero angle difference is always allowed!
        # More precisely minangle_rad <= 0 and maxaxangle_rad >= 0
        if branch.maxangle_rad < 0 or branch.minangle_rad > 0:
            raise ValueError(
                f"Broken assumption 1: branch j {j} f {branch.f} t {branch.t} "
                f"minanglerad {branch.minangle_rad} maxanglerad {branch.maxangle_rad}."
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
            lb=lbound, ub=ubound, name=f"c_{j}_{branch.f}_{branch.t}"
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
            lb=lbound, ub=ubound, name=f"s_{j}_{branch.f}_{branch.t}"
        )

    Pvar_f = {}
    Qvar_f = {}
    Pvar_t = {}
    Qvar_t = {}

    twinPvar_f = {}
    twinPvar_t = {}
    twinQvar_f = {}
    twinQvar_t = {}

    for j, branch in branches.items():
        busf = buses[IDtoCountmap[branch.f]]
        bust = buses[IDtoCountmap[branch.t]]

        if branch.constrainedflow:
            bound = branch.limit
        else:
            # Generous: assumes line charging up to 100%.
            # However it still amounts to an assumption.
            bound = 2 * (abs(alldata["summaxgenP"]) + abs(alldata["summaxgenQ"]))

        Pvar_f[branch] = model.addVar(
            lb=-bound, ub=bound, name=f"P_{j}_{branch.f}_{branch.t}"
        )

        Pvar_t[branch] = model.addVar(
            lb=-bound, ub=bound, name=f"P_{j}_{branch.t}_{branch.f}"
        )

        Qvar_f[branch] = model.addVar(
            lb=-bound, ub=bound, name=f"Q_{j}_{branch.f}_{branch.t}"
        )

        Qvar_t[branch] = model.addVar(
            lb=-bound, ub=bound, name=f"Q_{j}_{branch.t}_{branch.f}"
        )

        if alldata["branchswitching_mip"]:
            twinPvar_f[branch] = model.addVar(
                lb=-bound, ub=bound, name=f"twinP_{j}_{branch.f}_{branch.t}"
            )

            twinPvar_t[branch] = model.addVar(
                lb=-bound, ub=bound, name=f"twinP_{j}_{branch.t}_{branch.f}"
            )

            twinQvar_f[branch] = model.addVar(
                lb=-bound, ub=bound, name=f"twinQ_{j}_{branch.f}_{branch.t}"
            )

            twinQvar_t[branch] = model.addVar(
                lb=-bound, ub=bound, name=f"twinQ_{j}_{branch.t}_{branch.f}"
            )

    zvar = {}
    if alldata["branchswitching_mip"] or alldata["branchswitching_comp"]:
        for j, branch in branches.items():
            zvar[branch] = model.addVar(
                vtype=GRB.BINARY, name=f"z_{j}_{branch.f}_{branch.t}"
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
    # Pinjvar is the variable modeling total active power injected
    # by bus j into the branches incident with j
    alldata["LP"]["Pinjvar"] = Pinjvar
    # Qinjvar is the variable modeling total reactive power injected
    # by bus j into the branches incident with j
    alldata["LP"]["Qinjvar"] = Qinjvar


def lpformulator_ac_create_polar_vars(alldata, model):
    """
    Creates and adds variables for polar AC formulation to a given Gurobi model

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param: Gurobi model to be constructed
    :type model: :class: `gurobipy.Model`
    """

    buses = alldata["buses"]
    IDtoCountmap = alldata["IDtoCountmap"]
    branches = alldata["branches"]
    fixtolerance = alldata["fixtolerance"]
    newvarcount = 0

    logger.info("  Creating variables for polar formulation.")

    vvar = {}
    thetavar = {}
    for j, bus in buses.items():
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

    for j, branch in branches.items():
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

    buses = alldata["buses"]

    evar = {}
    fvar = {}

    for j, bus in buses.items():
        if alldata["usemaxdispersion"]:
            lbound = bus.Vmin * math.cos(alldata["maxdispersion_rad"])
            ubound = bus.Vmax
        else:
            ubound = bus.Vmax
            lbound = -ubound

        evar[bus] = model.addVar(lb=lbound, ub=ubound, name=f"e_{bus.nodeID}")

        if alldata["usemaxdispersion"]:
            lbound = 0
            ubound = bus.Vmax * math.sin(alldata["maxdispersion_rad"])
        elif j == alldata["refbus"]:
            ubound = lbound = 0

        fvar[bus] = model.addVar(lb=lbound, ub=ubound, name=f"f_{bus.nodeID}")

    # Save e, f variables data
    alldata["LP"]["evar"] = evar
    alldata["LP"]["fvar"] = fvar  # e^2 + f^2 = voltage magnitude


def lpformulator_ac_create_constraints(alldata, model):
    """
    Creates and adds constraints for AC formulation to a given Gurobi model

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param: Gurobi model to be constructed
    :type model: :class: `gurobipy.Model`
    """

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
    for j, branch in branches.items():
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
    for j, branch in branches.items():
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
    for j, bus in buses.items():
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

    for j, bus in buses.items():
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
    for j, bus in buses.items():
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
    for j, branch in branches.items():
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
        for j, branch in branches.items():
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
        for j, branch in branches.items():
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
            for j, branch in branches.items():
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
        for j, branch in branches.items():
            expr.add(zvar[branch])
        model.addConstr(expr >= N, name="sumz_lower_heuristic_bound")


def lpformulator_ac_add_polarconstraints(alldata, model):
    """
    Creates and adds constraints for polar represenation of AC formulation
    to a given Gurobi model

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param model: Gurobi model to be constructed
    :type model: :class: `gurobipy.Model`
    """

    buses = alldata["buses"]
    branches = alldata["branches"]
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

    for j, bus in buses.items():
        model.addConstr(cvar[bus] == vvar[bus] ** 2, name="cffdef_%d" % j)
        count += 1

    for j, branch in branches.items():
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
    Creates and adds nonconvex bilinear e, f constraints for AC formulation to
    a given Gurobi model

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param model: Gurobi model to be constructed
    :type model: :class: `gurobipy.Model`
    """

    buses = alldata["buses"]
    branches = alldata["branches"]
    evar = alldata["LP"]["evar"]
    fvar = alldata["LP"]["fvar"]  #
    # cvar[km] = (voltage mag at k) * (voltage mag at m) * cos(thetak - thetam)
    # Compute voltage magnitudes for k and m
    # Take arccos to get thetak - thetam = some value
    # Set arbitrarily theta1 = 0 and solve all others
    cvar = alldata["LP"]["cvar"]  # Square of the voltage magnitude e^2 + f^2
    svar = alldata["LP"]["svar"]
    IDtoCountmap = alldata["IDtoCountmap"]
    logger.info("  Adding nonconvex e, f, constraints.")

    for j, bus in buses.items():
        model.addConstr(
            -cvar[bus] + evar[bus] * evar[bus] + fvar[bus] * fvar[bus] == 0,
            name=f"cbusdef_{j}_{bus.nodeID}",
        )

    for j, branch in branches.items():
        if branch.status:
            busf = buses[IDtoCountmap[branch.f]]
            bust = buses[IDtoCountmap[branch.t]]

            model.addConstr(
                -cvar[branch] + evar[busf] * evar[bust] + fvar[busf] * fvar[bust] == 0,
                name=f"cdef_{j}_{branch.f}_{branch.t}",
            )
            model.addConstr(
                -svar[branch] - evar[busf] * fvar[bust] + fvar[busf] * evar[bust] == 0,
                name=f"sdef_{j}_{branch.f}_{branch.t}",
            )


def computebalbounds(alldata, bus):
    """
    Computes active and reactive max and min bus flow balance values

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
