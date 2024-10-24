import logging
import math

import gurobipy as gp

from gurobi_optimods.opf.grbformulator import OpfType, turn_solution_into_result_dict
from gurobi_optimods.opf.grbformulator_ac import lpformulator_ac_body

logger = logging.getLogger(__name__)


def compute_violations_from_voltages(env, alldata):
    """
    Calls the corresponding violation checker function

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param model: Gurobi model to be constructed
    :type model: :class: `gurobipy.Model`
    :param opftype: Type of OPF formulation
    :type opftype: :enum: `OpfType`

    :return: Dictionary holding case data following the MATPOWER notation with additional
             violations fields
    :rtype: dict
    """

    logger.info("Computing violations from given voltage inputs.")
    violations = None
    # Create model
    with gp.Model("AC_Violations_Model", env=env) as model:
        # Add model variables and constraints
        lpformulator_ac_body(alldata, model)
        # Compute violations
        lpformulator_ac_strictchecker(alldata, model)
        violations = turn_solution_into_result_dict(
            alldata, model, OpfType.AC, "violation"
        )

    return violations


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
    buses = alldata["buses"]
    IDtoCountmap = alldata["IDtoCountmap"]
    branches = alldata["branches"]
    cvar = alldata["LP"]["cvar"]
    svar = alldata["LP"]["svar"]

    for j, bus in buses.items():
        bus.inpute = bus.inputV * math.cos(bus.inputA_rad)
        bus.inputf = bus.inputV * math.sin(bus.inputA_rad)

    # Compute xbuffer
    if alldata["use_ef"]:
        evar = alldata["LP"]["evar"]
        fvar = alldata["LP"]["fvar"]
        for j, bus in buses.items():
            xbuffer[evar[bus]] = bus.inpute
            xbuffer[fvar[bus]] = bus.inputf
        logger.info("Derived e, f values.")

    if alldata["dopolar"] is False:
        for j, bus in buses.items():
            xbuffer[cvar[bus]] = bus.inpute * bus.inpute + bus.inputf * bus.inputf

        for j, branch in branches.items():
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

        for j, bus in buses.items():
            xbuffer[cvar[bus]] = bus.inputV * bus.inputV
            xbuffer[vvar[bus]] = bus.inputV
            xbuffer[thetavar[bus]] = bus.inputA_rad

        for j, branch in branches.items():
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

    # Recomputing the active and reactive power flows using the voltage variables
    for j, branch in branches.items():
        for item in [
            (branch.Pdeffconstr, Pvar_f[branch]),
            (branch.Pdeftconstr, Pvar_t[branch]),
            (branch.Qdeffconstr, Qvar_f[branch]),
            (branch.Qdeftconstr, Qvar_t[branch]),
        ]:
            constr = item[0]
            var = item[1]
            leadcoeff = 0
            # If the constraint is linear, e.g., Pvar_f = f(cosvar,sinvar), f linear
            if type(constr) is gp.Constr:
                row = model.getRow(constr)
                sum = -constr.RHS
                for i in range(row.size()):
                    v = row.getVar(i)
                    coeff = row.getCoeff(i)
                    if v.Varname != var.Varname:
                        sum += coeff * xbuffer[v]
                    else:
                        leadcoeff = coeff
            else:  # Rectangular formulation with Pvar_f = f(evar,fvar), f quadratic
                row = model.getQCRow(constr)
                sum = -constr.QCRHS
                for i in range(row.size()):
                    v1 = row.getVar1(i)
                    v2 = row.getVar2(i)
                    coeff = row.getCoeff(i)
                    sum += coeff * xbuffer[v1] * xbuffer[v2]
                lterms = row.getLinExpr()
                for i in range(lterms.size()):
                    v = lterms.getVar(i)
                    coeff = lterms.getCoeff(i)
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
    for j, bus in buses.items():
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
    branches = alldata["branches"]
    gens = alldata["gens"]

    max_violation_string = None
    max_violation_value = 0

    alldata["violation"] = {}  # Dictionary to keep track of violations
    Vmagviol = alldata["violation"]["Vmagviol"] = {}
    IPviol = alldata["violation"]["IPviol"] = {}
    IQviol = alldata["violation"]["IQviol"] = {}
    alldata["violation"]["branchlimit"] = {}
    for j, bus in buses.items():
        alldata["violation"][bus] = {}

    if alldata["use_ef"]:
        evar = alldata["LP"]["evar"]
        fvar = alldata["LP"]["fvar"]

    xbuffer = alldata["LP"]["xbuffer"]

    maxlbviol = maxubviol = 0
    badlbvar = None
    badubvar = None

    logger.info("Direct bus magnitude bounds check. Warning issued if violated.")

    for j, bus in buses.items():
        if bus.inputV > bus.Vmax:
            logger.warning(
                f">>> Warning: bus # {j} has input voltage {bus.inputV} "
                f"which is larger than Vmax {bus.Vmax}."
            )
            thisviol = bus.inputV - bus.Vmax
            if thisviol > max_violation_value:
                max_violation_string = "bus_" + str(bus.nodeID) + "_Vmax"
                max_violation_value = thisviol
        if bus.inputV < bus.Vmin:
            logger.warning(
                f">>> Warning: bus # {j} has input voltage {bus.inputV} "
                f"which is smaller than Vmin {bus.Vmin}."
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

    for j, branch in branches.items():
        fromvalue = math.sqrt(
            xbuffer[Pvar_f[branch]] * xbuffer[Pvar_f[branch]]
            + xbuffer[Qvar_f[branch]] * xbuffer[Qvar_f[branch]]
        )
        fromviol = max(fromvalue - branch.limit, 0)
        if fromvalue > branch.limit:
            logger.warning(
                f">>> Warning: branch # {j} has 'from' flow magnitude {fromvalue} "
                f"which is larger than limit {branch.limit}."
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
                f">>> Warning: branch # {j} has 'to' flow magnitude {tovalue} "
                f"which is larger than limit {branch.limit}."
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
        for j, bus in buses.items():
            for var in [(evar[bus], bus.inpute), (fvar[bus], bus.inputf)]:
                output_string = f"auxiliary variable {var[0].VarName} used for "
                f"rectangular formulation defined for bus ID {bus.nodeID}"
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
        for j, bus in buses.items():
            output_string = (
                f"variable {vvar[bus].VarName} defining voltage "
                f"magnitude of bus ID {bus.nodeID}"
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

    for j, branch in branches.items():
        for var in [Pvar_f[branch], Pvar_t[branch], Qvar_f[branch], Qvar_t[branch]]:
            power = (
                "real power injected"
                if "P" in var.VarName
                else "reactive power injected"
            )
            direction = "from" if "_f" in var.VarName else "to"
            output_string = (
                f"variable {var.VarName} defining {power} into '{direction}' "
                f"end of branch ({branch.f},{branch.t})"
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
    for j, bus in buses.items():
        varinj = Pinjvar[bus]
        output_string = (
            f"variable {varinj.VarName} defining total real power injected "
            f"by bus ID {bus.nodeID} into incident branches"
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

        # Construct min/max injection at the bus by looking at available
        # generators and load
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
    for j, bus in buses.items():
        varinj = Qinjvar[bus]
        varf = Qvar_f[branch]
        output_string = (
            f"variable {varinj.VarName} defining total reactive power "
            f"injected by bus ID {bus.nodeID} into incident branches"
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

        # but also, construct min/max injection at the bus by looking
        # at available generators and load
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
        f"Max overall violation {max_violation_value:.4e}, key: {max_violation_string}."
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
            f"LB violation for {output_string} LB {lb:<14.4e}  "
            f"x {value:<14.4e}  UB {ub:<14.4e}"
        )

    if lbviol > maxlbviol:
        maxlbviol = lb - value
        badlbvar = grbvariable

    if ubviol > 0:
        logger.info(
            f"UB violation for {output_string} LB {lb:<14.4e}  "
            f"x {value:<14.4e}  UB {ub:<14.4e}"
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
