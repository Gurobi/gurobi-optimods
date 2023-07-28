import logging

import gurobipy as gp
from gurobipy import GRB

from gurobi_optimods.opf.grbformulator_ac import lpformulator_ac_create_efvars
from gurobi_optimods.opf.grbformulator_common import set_gencost_objective

logger = logging.getLogger(__name__)


def lpformulator_iv_body(alldata, model):
    """
    Adds variables and constraints for IV formulation to a given Gurobi model

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param model: Gurobi model to be constructed
    :type model: :class: `gurobipy.Model`
    """

    # Assert compatible settings: Avoiding these asserts when called from the
    # public API should be handled correctly by build_internal_settings(...).
    # 1. Polar form incompatible with jabr
    if alldata["dopolar"]:
        assert alldata["skipjabr"]
    # 2. IV formulation currently does not support branch switching.
    assert not alldata["branchswitching_mip"]
    assert not alldata["branchswitching_comp"]
    # 3. IV formulation requires ef constraints.
    assert alldata["use_ef"]

    lpformulator_iv_create_vars(alldata, model)
    set_gencost_objective(alldata, model)
    lpformulator_iv_create_constraints(alldata, model)


def lpformulator_iv_create_vars(alldata, model):
    """
    Creates and adds variables for IV formulation to a given Gurobi model

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param model: Gurobi model to be constructed
    :type model: :class: `gurobipy.Model`
    """

    logger.info("Creating variables.")

    # fixtolerance = alldata["fixtolerance"] # Currently not used

    numbuses = alldata["numbuses"]
    buses = alldata["buses"]
    IDtoCountmap = alldata["IDtoCountmap"]
    gens = alldata["gens"]

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

        varcount += 1

        #at this point, minprod is the square of bus min voltage

    logger.info('Added %d cff variables\n'%(varcount))
    """

    # TODO shouldn't this check whether use_ef is True? Or is it always needed
    # (even with polar)?
    lpformulator_ac_create_efvars(alldata, model)

    # Next, generator variables
    GenPvar = {}
    GenQvar = {}

    for j in range(1, numbuses + 1):
        bus = buses[j]
        # First, injection variables
        bus.Vmax * bus.Vmax
        bus.Vmin * bus.Vmin

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
            buses[count_of_f].Vmax * buses[count_of_t].Vmax
            buses[count_of_f].Vmin * buses[count_of_t].Vmin
            if branch.constrainedflow:
                pass
            else:
                2 * (
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
            irvar_t[branch] = model.addVar(
                obj=0.0,
                lb=lbound,
                ub=ubound,
                name="ir_t_%d_%d_%d" % (j, bust.nodeID, busf.nodeID),
            )

            ijvar_f[branch] = model.addVar(
                obj=0.0,
                lb=lbound,
                ub=ubound,
                name="ij_f_%d_%d_%d" % (j, busf.nodeID, bust.nodeID),
            )
            ijvar_t[branch] = model.addVar(
                obj=0.0,
                lb=lbound,
                ub=ubound,
                name="ij_t_%d_%d_%d" % (j, bust.nodeID, busf.nodeID),
            )

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


def lpformulator_iv_create_constraints(alldata, model):
    """
    Creates and adds constraints for IV formulation to a given Gurobi model

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param model: Gurobi model to be constructed
    :type model: :class: `gurobipy.Model`
    """

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

    logger.info("Creating constraints.")

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

        logger.info(f"    {count} branch current definitions added.")

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
        logger.info(f"    {count} branch power flow definitions added.")
    else:  # Aggressive formulation
        logger.info(
            "Aggressive formulation, so will define power flows as functions of voltages."
        )

        logger.info("  Adding aggressive branch power definitions.")
        count = 0
        for j in range(1, 1 + numbranches):
            branch = branches[j]

            if not branch.status:
                continue

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
            qexpr.add(gkm * (evar[busf] * evar[bust] + fvar[busf] * fvar[bust]))
            qexpr.add(bkm * (-evar[busf] * fvar[bust] + fvar[busf] * evar[bust]))
            model.addConstr(qexpr == Pvar_f[branch], name=branch.Pfcname)

            # And imag power from k to m = -ek*[ bkk*ek + gkk*fk + bkm*em + gkm*fm ]
            #                            +  fk*[ gkk*ek - bkk*fk + gkm*em - bkm*fm ] =
            #
            #                 -bkk*(ek^2 + fk^2) - bkm*( ek*em + fk*fm) + gkm*(-ek*fm + em*fk)

            qexpr = -bkk * (evar[busf] * evar[busf] + fvar[busf] * fvar[busf])
            qexpr.add(-bkm * (evar[busf] * evar[bust] + fvar[busf] * fvar[bust]))
            qexpr.add(gkm * (-evar[busf] * fvar[bust] + fvar[busf] * evar[bust]))
            model.addConstr(qexpr == Qvar_f[branch], name=branch.Qfcname)

            # Now, the reversals

            gmm = branch.Gtt
            bmm = branch.Btt
            gmk = branch.Gtf
            bmk = branch.Btf

            qexpr = gmm * (evar[bust] * evar[bust] + fvar[bust] * fvar[bust])
            qexpr.add(gmk * (evar[bust] * evar[busf] + fvar[bust] * fvar[busf]))
            qexpr.add(bmk * (-evar[bust] * fvar[busf] + fvar[bust] * evar[busf]))
            model.addConstr(qexpr == Pvar_t[branch], name=branch.Ptcname)

            qexpr = -bmm * (evar[bust] * evar[bust] + fvar[bust] * fvar[bust])
            qexpr.add(-bmk * (evar[bust] * evar[busf] + fvar[bust] * fvar[busf]))
            qexpr.add(gmk * (-evar[bust] * fvar[busf] + fvar[bust] * evar[busf]))
            model.addConstr(qexpr == Qvar_t[branch], name=branch.Qtcname)

            count += 4

        logger.info(f"    {count} branch power flow definitions added.")

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

    logger.info(f"    {count} active bus power balance constraints added.")

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

    logger.info(f"    {count} reactive bus power balance constraints added.")

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
        model.addConstr(
            evar[bus] ** 2 + fvar[bus] ** 2 >= bus.Vmin * bus.Vmin, name="Vmin_%d" % j
        )
        count += 2

    logger.info(f"    {count} bus voltage limit constraints added.")

    # Branch limits.
    logger.info("  Adding branch limits.")
    count = 0
    for j in range(1, 1 + numbranches):
        branch = branches[j]

        if branch.status and branch.unboundedlimit is False:
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

    logger.info(f"    {count} branch limits added.")

    # Active loss inequalities.
    if alldata["useactivelossineqs"] is True:
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

        logger.info(f"    {count} active loss inequalities added.")
