import math
import time

import gurobipy as gp
from gurobipy import GRB
from myutils import break_exit
from grbfile import grbreadvoltsfile


def lpformulator_ac(alldata):
    """Formulate ACOPF model and solve it"""

    log = alldata['log']
    log.joint("\nAC formulation\n")

    starttime = time.time()

    # Handle special settings
    lpformulator_setup(alldata)

    # Create model
    with gp.Env() as env, gp.Model('grbacs', env=env) as model:
        # Add model variables and constraints
        lpformulator_ac_body(alldata, model)

        model.update() # Update to get correct model stats
        log.joint("Constructed ACOPF model with %d variables"%model.NumVars)
        log.joint(" and %d constraints\n\n"%model.NumConstrs)

        model.write(alldata['lpfilename'])  # FIXME remove.  Jarek: I am using this for debugging, for now
        log.joint('wrote LP to ' + alldata['lpfilename'] + '\n')
        break_exit('wrote lp') # 

        # Specific settings for better convergence
        model.params.NonConvex      = 2
        #model.params.DualReductions = 0

        #model.setParam(GRB.Param.MIPGap, 1.0e-10)
        #model.setParam(GRB.Param.FeasibilityTol, 1.0e-8)
        model.Params.MIPGap         = 1.0e-3
        model.Params.OptimalityTol  = 1.0e-3
        model.Params.FeasibilityTol = 1.0e-6 # 1.0e-8

        feastol = model.Params.FeasibilityTol
        opttol  = model.Params.OptimalityTol
        mipgap  = model.Params.MIPGap

        # Optimize
        model.optimize()

        # Check model status and re-optimize or try computing an IIS if necessary
        if model.status == GRB.INF_OR_UNBD:
            log.joint("\nModel Status: infeasible or unbounded\n")
            log.joint("Re-optimizing with DualReductions turned off\n\n")
            model.Params.DualReductions = 0
            model.optimize()

        if model.status == GRB.INFEASIBLE:
            log.joint("\nModel Status: infeasible")
            log.joint("Computing IIS...\n\n")
            model.computeIIS()
            log.joint("\nIIS computed, writing IIS to file acopfmodel.ilp\n\n")
            model.write("acopfmodel.ilp")

        elif model.status == GRB.UNBOUNDED:
            log.joint("\nModel Status: unbounded\n\n")

        elif model.status == GRB.INTERRUPTED:
            log.joint("\nModel Status: interrupted\n\n")

        elif model.status == GRB.OPTIMAL:
            log.joint("\nModel Status: optimal\n\n")

        # Only print objective value and solution quality if at least
        # one feasible point is available
        if model.SolCount > 0:
            log.joint("Objective value = %g"%model.objVal)
            model.printQuality()
            # FIXME yes, to be done
            # Here we should gather optimal solution values and gather them
            # in a standardized format which ultimately will be returned to the user
            # 

    endtime = time.time()
    log.joint("Overall time taken (model construction + optimization): %f s\n"%(endtime-starttime))

    '''
    buses        = alldata['buses']
    branches     = alldata['branches']
    gens         = alldata['gens']
    IDtoCountmap = alldata['IDtoCountmap']
    '''

def lpformulator_setup(alldata):
    """Helper function to handle specific settings"""

    alldata['maxdispersion_rad'] = (math.pi/180.)*alldata['maxdispersion_deg']

    if alldata['dopolar']:
        log = alldata['log']
        log.joint("  Polar formulation, shutting down incompatible options:\n")
        alldata['use_ef']               = False
        alldata['useconvexformulation'] = False
        alldata['skipjabr']             = True
        log.joint("    use_ef, useconvexformulation, jabr\n")

    if alldata['voltsfilename'] != 'NONE':
        grbreadvoltsfile(alldata)

def lpformulator_ac_body(alldata, model):
    """Helper function for adding variables and constraints to the model"""

    # Create model variables
    lpformulator_create_ac_vars(alldata, model)
    # Create model constraints
    lpformulator_create_ac_constraints(alldata, model)

def lpformulator_create_ac_vars(alldata, model):
    """Create model variables for ACOPF"""

    log = alldata['log']
    log.joint('Creating variables\n')

    fixtolerance = 1e-05

    if alldata['fixtolerance'] > 0:
        fixtolerance = alldata['fixtolerance']

    numbuses     = alldata['numbuses']
    buses        = alldata['buses']
    IDtoCountmap = alldata['IDtoCountmap']

    # Bus related variables
    cvar     = {}
    svar     = {}
    Pinjvar  = {}
    Qinjvar  = {}
    GenPvar  = {}
    GenQvar  = {}
    gens     = alldata['gens']
    varcount = 0



    for j in range(1,numbuses+1):
        bus = buses[j]

        # First, injection variables
        maxprod = bus.Vmax*bus.Vmax
        minprod = bus.Vmin*bus.Vmin
        ubound  = maxprod
        lbound  = minprod

        if alldata['FIXCS'] and bus.inputvoltage:
            lbound = bus.inputV*bus.inputV - fixtolerance
            ubound = bus.inputV*bus.inputV + fixtolerance

        cvar[bus] = model.addVar(obj = 0.0, lb = lbound, ub = ubound,
                                 name = "c_%d_%d"%(bus.nodeID, bus.nodeID))

        bus.cffvarind = varcount
        varcount += 1

        # csdefslacks to be done
        Plbound = Qlbound = -GRB.INFINITY
        Pubound = Qubound = GRB.INFINITY
        Pubound, Plbound, Qubound, Qlbound = computebalbounds(log, alldata, bus)

        Pinjvar[bus] = model.addVar(obj = 0.0, lb = Plbound, ub = Pubound,
                                    name = "IP_%d"%bus.nodeID)
        bus.Pinjvarind = varcount
        varcount += 1
        Qinjvar[bus] = model.addVar(obj = 0.0, lb = Qlbound, ub = Qubound,
                                    name = "IQ_%d"%bus.nodeID)
        bus.Qinjvarind = varcount
        varcount += 1

        # Next, generator variables
        for genid in bus.genidsbycount:
            gen   = gens[genid]
            lower = gen.Pmin*gen.status
            upper = gen.Pmax*gen.status
            #if bus.nodetype == 3:
            #  upper = GRB.INFINITY
            #  lower = -GRB.INFINITY  #ignoring slack bus
            GenPvar[gen] = model.addVar(obj = 0.0, lb = lower, ub = upper,
                                        name = "GP_%d_%d"%(gen.count,gen.nodeID))
            gen.Pvarind = varcount
            varcount += 1

            lower = gen.Qmin*gen.status
            upper = gen.Qmax*gen.status
            if bus.nodetype == 3:
                lower = -GRB.INFINITY
                upper = GRB.INFINITY

            GenQvar[gen] = model.addVar(obj = 0.0, lb = lower, ub = upper,
                                        name = "GQ_%d_%d"%(gen.count,gen.nodeID))
            gen.Qvarind = varcount
            varcount += 1


    # Branch related variables
    branches    = alldata['branches']
    numbranches = alldata['numbranches']
    for j in range(1,1+numbranches):
        branch     = branches[j]
        f          = branch.f
        t          = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]
        busf       = buses[count_of_f]
        bust       = buses[count_of_t]
        maxprod    = buses[count_of_f].Vmax*buses[count_of_t].Vmax
        minprod    = buses[count_of_f].Vmin*buses[count_of_t].Vmin

        # Assumption 1.  zero angle difference is always allowed! More precisely minangle_rad <= 0 and maxaxangle_rad >= 0
        if branch.maxangle_rad < 0 or branch.minangle_rad > 0:
            log.joint(' --- Broken assumption 1: branch j %d f %d t %d minanglerad %f maxanglerad %f\n'%(j, f, t, branch.minangle_rad, branch.maxangle_rad))
            log.raise_exception("phase angle assumption 1\n")

        ubound = ubasic = maxprod
        lbound = lbasic = -maxprod
        maxanglerad = branch.maxangle_rad
        minanglerad = branch.minangle_rad

        # Cosine
        if maxanglerad <= 0.5*math.pi:
            # In this case minangle <= 0
            if minanglerad >= -0.5*math.pi:
                lbound = minprod*min(math.cos(maxanglerad), math.cos(minanglerad))
            elif minanglerad >= -math.pi:
                lbound = maxprod*math.cos(minangle_rad)  # Which is negative
            elif minanglerad >= -1.5*math.pi:
                lbound = -maxprod
            else:
                lbound = -maxprod

        elif maxanglerad <= math.pi:
            if minanglerad >= -0.5*math.pi:
                lbound = maxprod*math.cos(maxanglerad)
            elif minanglerad >= -math.pi:
                lbound = maxprod*min(math.cos(maxanglerad), math.cos(minanglerad))
            elif minanglerad >= -1.5*math.pi:
                lbound = -maxprod
            else:
                lbound = -maxprod

        elif maxanglerad <= 1.5*math.pi:
            lbound = -maxprod

        elif maxanglerad <= 2*math.pi:
            lbound = -maxprod

        else:
            ubound = maxprod
            lbound = -maxprod

        if branch.inputcs:
            ubound = branch.inputc + fixtolerance
            lbound = branch.inputc - fixtolerance

        cvar[branch] = model.addVar(obj = 0.0, lb = lbound, ub = ubound,
                                    name = "c_%d_%d_%d"%(j, busf.nodeID, bust.nodeID))
        branch.cftvarind = varcount
        varcount += 1

        # Sine
        if maxanglerad <= math.pi/2:
            ubound = maxprod*sin(maxanglerad)

            if  minanglerad >= -0.5*math.pi:
                lbound = maxprod*sin(minanglerad)
            elif  minanglerad >= -math.pi:
                lbound = -maxprod
            elif  minanglerad >= -1.5*math.pi:
                ubound = maxprod*max( sin(maxanglerad), sin(minanglerad))
                lbound = -maxprod
            else:
                ubound = maxprod
                lbound = -maxprod

        elif maxanglerad <= math.pi:
            ubound = maxprod

            if minanglerad >= -0.5*math.pi:
                lbound = maxprod*sin(minanglerad)
            elif minanglerad >= -math.pi:
                lbound = -maxprod
            elif minanglerad >= -1.5*math.pi:
                lbound = -maxprod
            else:
                lbound = -maxprod

        elif maxanglerad <= 1.5*math.pi:
            ubound = maxprod

            if minanglerad >= -0.5*math.pi:
                lbound = maxprod*min(sin(maxanglerad), sin(minanglerad))
            elif minanglerad >= -math.pi:
                lbound = -maxprod
            elif minanglerad >= -1.5*math.pi:
                lbound = -maxprod
            else:
                lbound = -maxprod
        else:
            ubound = maxprod
            lbound = -maxprod

        svar[branch] = model.addVar(obj = 0.0, lb = lbound, ub = ubound,
                                    name = "s_%d_%d_%d"%(j, busf.nodeID, bust.nodeID))
        branch.sftvarind = varcount
        varcount += 1

    Pvar_f = {}
    Qvar_f = {}
    Pvar_t = {}
    Qvar_t = {}

    for j in range(1,1+numbranches):
        branch     = branches[j]
        f          = branch.f
        t          = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]
        busf       = buses[count_of_f]
        bust       = buses[count_of_t]
        ubound     = branch.limit
        lbound     = -ubound

        Pvar_f[branch] = model.addVar(obj = 0.0, lb = lbound, ub = ubound,
                                      name = "P_%d_%d_%d"%(j, busf.nodeID, bust.nodeID))
        branch.Pftvarind = varcount
        varcount += 1

        Pvar_t[branch] = model.addVar(obj = 0.0, lb = lbound, ub = ubound,
                                      name = "P_%d_%d_%d"%(j, bust.nodeID, busf.nodeID))
        branch.Ptfvarind = varcount
        varcount += 1

        Qvar_f[branch] = model.addVar(obj = 0.0, lb = lbound, ub = ubound,
                                      name = "Q_%d_%d_%d"%(j, busf.nodeID, bust.nodeID))
        branch.Qftvarind = varcount
        varcount += 1

        Qvar_t[branch] = model.addVar(obj = 0.0, lb = lbound, ub = ubound,
                                      name = "Q_%d_%d_%d"%(j, bust.nodeID, busf.nodeID))
        branch.Qtfvarind = varcount
        varcount += 1


    zvar = {}


    if alldata['branchswitching_mip'] or alldata['branchswitching_comp']:
        log.joint('Adding branch switching variables.\n')
        for j in range(1,1+numbranches):
            branch     = branches[j]
            f          = branch.f
            t          = branch.t
            lbound = 0
            ubound = 1
            zvar[branch] = model.addVar(lb = lbound, ub = ubound, obj = 0.0, vtype = GRB.INTEGER,
                                      name = "z_%d_%d_%d"%(j, f, t))
            branch.switchvarind = varcount
            varcount += 1
    alldata['LP']['zvar'] = zvar

    lincostvar = model.addVar(obj = 1.0, lb = -GRB.INFINITY, ub = GRB.INFINITY,
                              name = "lincost")
    alldata['LP']['lincostvar']    = lincostvar
    alldata['LP']['lincostvarind'] = varcount
    varcount += 1

    if alldata['usequadcostvar']:
        quadcostvar = model.addVar(obj = 1.0, lb = 0, ub = GRB.INFINITY,
                                   name = "quadcost")
        alldata['LP']['quadcostvar']    = quadcostvar
        alldata['LP']['quadcostvarind'] = varcount
        varcount += 1

    constobjval = 0
    for gen in gens.values():
        if gen.status > 0:
            constobjval += gen.costvector[gen.costdegree]

    constvar = model.addVar(obj = constobjval, lb = 1.0, ub = 1.0,
                            name = "constant")
    alldata['LP']['constvar'] = constvar
    varcount += 1

    # Powerflow variables
    if alldata['use_ef'] and alldata['useconvexformulation'] == False:
        lpformulator_create_ac_efvars(alldata, model, varcount)

    if alldata['dopolar']:
        lpformulator_create_ac_polar_vars(alldata, model, varcount)

    # Save variable data
    alldata['LP']['cvar']   = cvar
    alldata['LP']['svar']   = svar
    alldata['LP']['Pvar_f'] = Pvar_f
    alldata['LP']['Pvar_t'] = Pvar_t
    alldata['LP']['Qvar_f'] = Qvar_f
    alldata['LP']['Qvar_t'] = Qvar_t
    alldata['LP']['GenPvar'] = GenPvar
    alldata['LP']['GenQvar'] = GenQvar
    alldata['LP']['Pinjvar'] = Pinjvar
    alldata['LP']['Qinjvar'] = Qinjvar

def lpformulator_create_ac_polar_vars(alldata, model, varcount):
    """Create polar variables for ACOPF"""

    vvar         = {}
    thetavar     = {}
    numbuses     = alldata['numbuses']
    buses        = alldata['buses']
    IDtoCountmap = alldata['IDtoCountmap']
    numbranches  = alldata['numbranches']
    branches     = alldata['branches']
    log          = alldata['log']
    newvarcount  = 0
    fixtolerance = alldata['fixtolerance']

    log.joint("  Creating variables for polar formulation\n")

    for j in range(1,numbuses+1):
        bus       = buses[j]
        ubound    = bus.Vmax
        lbound    = bus.Vmin

        if bus.inputvoltage:
            candidatelbound = bus.inputV - fixtolerance
            candidateubound = bus.inputV + fixtolerance

            lbound = max(lbound, candidatelbound)
            ubound = min(ubound, candidateubound)
            #print('------>',j, bus.nodeID, bus.inputV, 'lbound',lbound, 'ubound',ubound)
        vvar[bus] = model.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "v_"+str(bus.nodeID))
        ubound = 2*math.pi
        lbound = -ubound

        if bus.inputvoltage:
            candidatelbound = bus.inputA_rad - fixtolerance
            candidateubound = bus.inputA_rad + fixtolerance
            lbound = max(lbound, candidatelbound)
            ubound = min(ubound, candidateubound)

    
        thetavar[bus]  = model.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "theta_"+str(bus.nodeID))
        newvarcount += 2


    cosvar     = {}
    sinvar     = {}
    thetaftvar = {}
    vfvtvar    = {}

    log.joint("    Assumption. Phase angle diffs between -pi and pi\n")

    for j in range(1,1+numbranches):
        branch     = branches[j]
        f          = branch.f
        t          = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]
        busf       = buses[count_of_f]
        bust       = buses[count_of_t]

        cosvar[branch] = model.addVar(obj = 0.0, lb = -1.0, ub = 1.0,
                                      name = "cos_%d_%d_%d"%(j, busf.nodeID, bust.nodeID))
        sinvar[branch] = model.addVar(obj = 0.0, lb = -1.0, ub = 1.0,
                                      name = "sin_%d_%d_%d"%(j, busf.nodeID, bust.nodeID))

        # Major assumption! Heuristic: angle diffs between -pi and pi.
        thetaftvar[branch] = model.addVar(obj = 0.0, lb = -math.pi, ub = math.pi,
                                          name = "thetaft_%d_%d_%d"%(j, busf.nodeID, bust.nodeID))

        ubound = busf.Vmax*bust.Vmax
        lbound = busf.Vmin*bust.Vmin

        vfvtvar[branch] = model.addVar(obj = 0.0, lb = lbound, ub = ubound,
                                       name = "vfvt_%d_%d_%d"%(j, busf.nodeID, bust.nodeID))
        newvarcount += 4

    # Save polar variables data
    alldata['LP']['vvar']       = vvar
    alldata['LP']['thetavar']   = thetavar
    alldata['LP']['cosvar']     = cosvar
    alldata['LP']['sinvar']     = sinvar
    alldata['LP']['thetaftvar'] = thetaftvar
    alldata['LP']['vfvtvar']    = vfvtvar

    log.joint("    Added %d new variables to handle polar formulation.\n"%newvarcount)
    break_exit('polarvars')

    varcount += newvarcount

def lpformulator_create_ac_efvars(alldata, model, varcount):
    """Create nonconvex e, f variables for ACOPF"""

    evar         = {}
    fvar         = {}
    numbuses     = alldata['numbuses']
    buses        = alldata['buses']
    IDtoCountmap = alldata['IDtoCountmap']
    efvarcount   = 0
    log          = alldata['log']
    fixtolerance = 1e-05
    if alldata['fixtolerance'] > 0:
        fixtolerance = alldata['fixtolerance']

    log.joint("  Creating e,f variables\n")

    for j in range(1,numbuses+1):
        bus    = buses[j]
        ubound = bus.Vmax
        lbound = -ubound

        if alldata['usemaxdispersion']:
            lbound = bus.Vmin*math.cos(alldata['maxdispersion_rad'])
            ubound = bus.Vmax

        evar[bus] = model.addVar(obj = 0.0, lb = lbound, ub = ubound,
                                 name = "e_%d"%bus.nodeID)
        bus.evarind = varcount + efvarcount
        efvarcount += 1

        if alldata['usemaxdispersion']:
            lbound = 0
            ubound = bus.Vmax*math.sin(alldata['maxdispersion_rad'])
        elif j == alldata['refbus']:
            ubound = lbound = 0

        fvar[bus] = model.addVar(obj = 0.0, lb = lbound, ub = ubound,
                                 name = "f_%d"%bus.nodeID)
        bus.fvarind = varcount + efvarcount
        efvarcount += 1

    # Save e, f variables data
    alldata['LP']['evar'] = evar
    alldata['LP']['fvar'] = fvar

    log.joint("  Added %d e, f variables\n"%efvarcount)

    varcount += efvarcount

def computebalbounds(log, alldata, bus):
    """Compute bal bounds"""
    # Computes active and reactive max and min bus flow balance values
    
    # First let's get max/min generations
    loud = 0 # See below

    baseMVA = alldata['baseMVA']
    gens    = alldata['gens']
    Pubound = Plbound = 0
    Qubound = Qlbound = 0

    for gencounter in bus.genidsbycount:
        if gens[gencounter].status:
            Pubound += gens[gencounter].Pmax
            Plbound += gens[gencounter].Pmin
            Qubound += gens[gencounter].Qmax
            Qlbound += gens[gencounter].Qmin

        if loud: # it is worth keeping because a user may want to debug the generator limits that they imposed, which could make a problem infeasible
            #log.joint(" Pubound for %d %f genc %d\n"%(bus.nodeID, Pubound, gencounter))
            #log.joint(" Plbound for %d %f genc %d\n"%(bus.nodeID, Plbound, gencounter))
            log.joint(" Qubound for %d %f genc %d\n"%(bus.nodeID, Qubound, gencounter))
            log.joint(" Qlbound for %d %f genc %d\n"%(bus.nodeID, Qlbound, gencounter))

    Pubound -= bus.Pd
    Plbound -= bus.Pd
    Qubound -= bus.Qd
    Qlbound -= bus.Qd

    if bus.nodetype == 4:
        Pubound = Plbound = Qubound = Qlbound = 0

    if loud: # same as above
        #log.joint(" Pubound for %d final %f\n"%(bus.nodeID, Pubound))
        #log.joint(" (Pd was %g)\n"%bus.Pd)
        #log.joint(" Plbound for %d final %f\n"%(bus.nodeID, Plbound))
        log.joint(" Qubound for %d final %f\n"%(bus.nodeID, Qubound))
        log.joint(" (Qd was %g)\n"%bus.Qd)
        log.joint(" Qlbound for %d final %f\n"%(bus.nodeID, Qlbound))
        break_exit(" ")

    return Pubound, Plbound, Qubound, Qlbound

def lpformulator_create_ac_constraints(alldata, model):
    """"Create constraint for ACOPF"""

    numbuses     = alldata['numbuses']
    buses        = alldata['buses']
    numbranches  = alldata['numbranches']
    branches     = alldata['branches']
    gens         = alldata['gens']
    IDtoCountmap = alldata['IDtoCountmap']
    cvar         = alldata['LP']['cvar']
    svar         = alldata['LP']['svar']
    Pvar_f       = alldata['LP']['Pvar_f']
    Pvar_t       = alldata['LP']['Pvar_t']
    Qvar_f       = alldata['LP']['Qvar_f']
    Qvar_t       = alldata['LP']['Qvar_t']
    Pinjvar      = alldata['LP']['Pinjvar']
    Qinjvar      = alldata['LP']['Qinjvar']
    GenPvar      = alldata['LP']['GenPvar']
    GenQvar      = alldata['LP']['GenQvar']
    lincostvar   = alldata['LP']['lincostvar']
    zvar         = alldata['LP']['zvar']
    log          = alldata['log']

    log.joint("Creating constraints\n")
    log.joint("  Adding cost definition constraint\n")

    coeff     = [gen.costvector[gen.costdegree - 1] for gen in gens.values()]
    variables = [GenPvar[gen] for gen in gens.values()]
    expr      = gp.LinExpr(coeff, variables)
    model.addConstr(expr == lincostvar, name = "lincostdef")

    numquadgens = 0
    for gen in gens.values():
        if gen.costdegree >= 2 and gen.costvector[0] > 0 and gen.status:
            numquadgens += 1

    log.joint("    Number of generators with quadratic cost coefficient: %d\n"%numquadgens)

    if numquadgens > 0:
        if alldata['usequadcostvar']:
            quadcostvar = alldata['LP']['quadcostvar']
            log.joint("    Adding quadcost definition constraint\n")
            qcost = gp.QuadExpr()
            for gen in gens.values():
                if gen.costdegree == 2 and gen.costvector[0] != 0:
                    qcost.add(gen.costvector[0]*GenPvar[gen]*GenPvar[gen])

            model.addConstr(qcost <= quadcostvar, name = "qcostdef")
        else:
            log.joint("    Adding quad cost to objective\n")
            model.update() # Necessary to flush changes in the objective function
            oldobj = model.getObjective()
            newobj = gp.QuadExpr(oldobj)
            for gen in gens.values():
                if gen.costdegree == 2 and gen.costvector[0] != 0:
                    newobj.add(gen.costvector[0]*GenPvar[gen]*GenPvar[gen])

            model.setObjective(newobj, GRB.MINIMIZE)

    # Active PF defs
    log.joint("  Adding active power flow definitions\n")
    count = 0
    for j in range(1,1+numbranches):
        branch     = branches[j]
        f          = branch.f
        t          = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]
        busf       = buses[count_of_f]
        bust       = buses[count_of_t]
        if branch.status:
            if alldata['substitute_nonconv'] == False or alldata['use_ef'] == False:
                # Gff cff + Gft cft + Bft sft
                expr = gp.LinExpr([branch.Gff, branch.Gft, branch.Bft],
                              [cvar[busf], cvar[branch], svar[branch]])
                model.addConstr(expr == Pvar_f[branch], name = "Pdef_%d_%d_%d"%(j, f, t))
                # Gtt ctt + Gtf cft + Btf stf = Gtt ctt + Gtf cft - Btf sft
                expr = gp.LinExpr([branch.Gtt, branch.Gtf, -branch.Btf], # minus because svarft = -svartf
                                  [cvar[bust], cvar[branch], svar[branch]])
                model.addConstr(expr == Pvar_t[branch], name = "Pdef_%d_%d_%d"%(j, t, f))
                count += 2
        else:
            model.addConstr(Pvar_f[branch] == 0, name = "Pdef_%d_%d_%d"%(j, f, t))
            model.addConstr(Pvar_t[branch] == 0, name = "Pdef_%d_%d_%d"%(j, t, f))            

    log.joint("    %d active power flow definitions added.\n"%count)

    # Reactive PF defs
    log.joint("  Adding reactive power flow definitions\n")
    count = 0
    for j in range(1,1+numbranches):
        branch     = branches[j]
        f          = branch.f
        t          = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]
        busf       = buses[count_of_f]
        bust       = buses[count_of_t]

        if branch.status:
            if alldata['substitute_nonconv'] == False or alldata['use_ef'] == False:
                # -Bff cff - Bft cft + Gft sft
                expr = gp.LinExpr([-branch.Bff, -branch.Bft, branch.Gft],
                              [cvar[busf], cvar[branch], svar[branch]])
                model.addConstr(expr == Qvar_f[branch], name = "Qdef_%d_%d_%d"%(j, f, t))
                # -Btt ctt - Btf cft + Gtf stf = -Btt ctt - Btf cft - Gtf sft 
                expr = gp.LinExpr([-branch.Btt, -branch.Btf, -branch.Gtf], # again, same minus
                              [cvar[bust], cvar[branch], svar[branch]])
                model.addConstr(expr == Qvar_t[branch], name = "Qdef_%d_%d_%d"%(j, t, f))
                count += 2
        else:
            model.addConstr(Qvar_f[branch] == 0, name = "Qdef_%d_%d_%d"%(j, f, t))
            model.addConstr(Qvar_t[branch] == 0, name = "Qdef_%d_%d_%d"%(j, t, f))            

    log.joint("    %d reactive power flow definitions added\n"%count)

    # Balance constraints
    log.joint("  Adding balance constraints\n")
    count = 0
    for j in range(1,1+numbuses):
        bus  = buses[j]
        expr = gp.LinExpr()
        qexpr = gp.QuadExpr()
        if bus.Gs != 0:
            expr.add(bus.Gs*cvar[bus])

        if alldata['branchswitching_comp'] == False:
            for branchid in bus.frombranchids.values():
                expr.add(Pvar_f[branches[branchid]])

            for branchid in bus.tobranchids.values():
                expr.add(Pvar_t[branches[branchid]])
        else:
            #product of binary and powerflow
            for branchid in bus.frombranchids.values():
                qexpr.add(zvar[branches[branchid]]*Pvar_f[branches[branchid]])

            for branchid in bus.tobranchids.values():
                qexpr.add(zvar[branches[branchid]]*Pvar_t[branches[branchid]])
            
        model.addConstr(expr + qexpr == Pinjvar[bus], name = "PBaldef%d_%d"%(j, bus.nodeID))
        
        count += 1

    for j in range(1,1+numbuses):
        bus  = buses[j]
        expr = gp.LinExpr()
        qexpr = gp.QuadExpr()
        
        if bus.Bs != 0:
            expr.add(-bus.Bs * cvar[bus])
        if alldata['branchswitching_comp'] == False:
            for branchid in bus.frombranchids.values():
                expr.add(Qvar_f[branches[branchid]])

            for branchid in bus.tobranchids.values():
                expr.add(Qvar_t[branches[branchid]])
        else:
            #product of binary and powerflow
            for branchid in bus.frombranchids.values():
                qexpr.add(zvar[branches[branchid]]*Qvar_f[branches[branchid]])

            for branchid in bus.tobranchids.values():
                qexpr.add(zvar[branches[branchid]]*Qvar_t[branches[branchid]])

        model.addConstr(expr + qexpr == Qinjvar[bus], name = "QBaldef%d_%d"%(j, bus.nodeID))
        count += 1

    log.joint("    %d balance constraints added\n"%count)

    # Injection defs
    log.joint("  Adding injection definition constraints\n")
    count = 0
    for j in range(1,1+numbuses):
        bus  = buses[j]
        expr = gp.LinExpr()

        if len(bus.genidsbycount) > 0:
            for genid in bus.genidsbycount:
                gen = gens[genid]
                expr.add(GenPvar[gen])

        model.addConstr(Pinjvar[bus] == expr - bus.Pd, name = "Bus_PInj_%d"%j)

        expr = gp.LinExpr()

        if len(bus.genidsbycount) > 0:
            for genid in bus.genidsbycount:
                gen = gens[genid]
                expr.add(GenQvar[gen])

        model.addConstr(Qinjvar[bus] == expr - bus.Qd, name = "Bus_QInj_%d"%j)
        count += 2

    log.joint("    %d injection definition constraints added\n"%count)

    # Branch limits
    log.joint("  Adding branch limits\n")
    count = 0
    for j in range(1,1+numbranches):
        branch = branches[j]

        if branch.status and branch.unboundedlimit == False:
            f          = branch.f
            t          = branch.t
            count_of_f = IDtoCountmap[f]
            count_of_t = IDtoCountmap[t]
            busf       = buses[count_of_f]
            bust       = buses[count_of_t]
            model.addConstr(Pvar_f[branch]*Pvar_f[branch] + Qvar_f[branch]*Qvar_f[branch] <= branch.limit**2,
                            name = "limit_f_%d_%d_%d"%(j, f, t))
            #themodel.cbLazy(Pvar_t[branch]*Pvar_t[branch] + Qvar_t[branch]*Qvar_t[branch] <= branch.limit**2)
            model.addConstr(Pvar_t[branch]*Pvar_t[branch] + Qvar_t[branch]*Qvar_t[branch] <= branch.limit**2,
                            name = "limit_t_%d_%d_%d"%(j, t, f))
            count += 2

    log.joint("    %d branch limits added\n"%count)

    # JABR
    if alldata['skipjabr'] == False:
        log.joint("  Adding Jabr constraints\n")
        count = 0
        for j in range(1,1+numbranches):
            branch = branches[j]
            if branch.status:
                f          = branch.f
                t          = branch.t
                count_of_f = IDtoCountmap[f]
                count_of_t = IDtoCountmap[t]
                busf       = buses[count_of_f]
                bust       = buses[count_of_t]
                model.addConstr(cvar[branch]*cvar[branch] + svar[branch]*svar[branch] <= cvar[busf]*cvar[bust],
                                name = "jabr_%d_%d_%d"%(j, f, t))
                count += 1

        log.joint("    %d Jabr constraints added\n"%count)
    else:
        log.joint("  Skipping Jabr inequalities\n")

    # Polar representation
    if alldata['dopolar']:
        lpformulator_ac_add_polarconstraints(alldata, model)

    # nonconvex e, f representation
    if alldata['use_ef'] and alldata['useconvexformulation'] == False:
        lpformulator_ac_add_nonconvexconstraints(alldata, model)

def lpformulator_ac_add_polarconstraints(alldata,model):
    """Create polar representation constraints for ACOPF"""

    buses        = alldata['buses']
    numbuses     = alldata['numbuses']
    branches     = alldata['branches']
    numbranches  = alldata['numbranches']
    vvar         = alldata['LP']['vvar']
    thetavar     = alldata['LP']['thetavar']
    thetaftvar   = alldata['LP']['thetaftvar']
    vfvtvar      = alldata['LP']['vfvtvar']
    cosvar       = alldata['LP']['cosvar']
    sinvar       = alldata['LP']['sinvar']
    IDtoCountmap = alldata['IDtoCountmap']
    cvar         = alldata['LP']['cvar']
    svar         = alldata['LP']['svar']
    log          = alldata['log']
    count        = 0

    log.joint("  Adding polar constraints\n")

    for j in range(1,1+numbuses):
        bus = buses[j]
        model.addConstr(cvar[bus] == vvar[bus]**2, name = "cffdef_%d"%j)
        count += 1

    for j in range(1,1+numbranches):
        branch     = branches[j]
        f          = branch.f
        t          = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]
        busf       = buses[count_of_f]
        bust       = buses[count_of_t]
        model.addConstr(thetaftvar[branch] == thetavar[busf] - thetavar[bust],
                        name = "thetaftdef_%d"%j)

        gc = model.addGenConstrCos(thetaftvar[branch], cosvar[branch],
                                   name = "cosdef_%d_%d_%d"%(j, f, t))
        gs = model.addGenConstrSin(thetaftvar[branch], sinvar[branch],
                                   name = "sindef_%d_%d_%d"%(j, f, t))
        model.addConstr(vfvtvar[branch] == vvar[busf]*vvar[bust],
                        name = "vfvtdef_%d_%d_%d"%(j, f, t))

        model.addConstr(cvar[branch] == vfvtvar[branch]*cosvar[branch],
                        name = "cftdef_%d"%j)
        model.addConstr(svar[branch] == vfvtvar[branch]*sinvar[branch],
                        name = "sftdef_%d"%j)
        count += 6 # Count general constraints as well

    log.joint("    %d polar constraints added\n"%count)

    break_exit('polarconstrs')

def lpformulator_ac_add_nonconvexconstraints(alldata,model):
    """Create nonconvex e, f constraints"""

    buses        = alldata['buses']
    numbuses     = alldata['numbuses']
    branches     = alldata['branches']
    numbranches  = alldata['numbranches']
    evar         = alldata['LP']['evar']
    fvar         = alldata['LP']['fvar']
    cvar         = alldata['LP']['cvar']
    svar         = alldata['LP']['svar']
    IDtoCountmap = alldata['IDtoCountmap']
    log          = alldata['log']
    count        = 0
    log.joint("  Adding nonconvex e, f, constraints\n")

    for j in range(1,1+numbuses):
        bus = buses[j]
        model.addConstr(-cvar[bus] + evar[bus]*evar[bus] + fvar[bus]*fvar[bus] == 0,
                        name = "cbusdef_%d_%d"%(j, bus.nodeID))
        count += 1

    for j in range(1,1+numbranches):
        branch = branches[j]

        if branch.status:
            f          = branch.f
            t          = branch.t
            count_of_f = IDtoCountmap[f]
            count_of_t = IDtoCountmap[t]
            busf       = buses[count_of_f]
            bust       = buses[count_of_t]

            model.addConstr(-cvar[branch] + evar[busf]*evar[bust] + fvar[busf]*fvar[bust] == 0,
                            name = "cdef_%d_%d_%d"%(j, f, t))
            model.addConstr(-svar[branch] - evar[busf]*fvar[bust] + fvar[busf]*evar[bust] == 0,
                            name = "sdef_%d_%d_%d"%(j, f, t))
            count += 2

    log.joint("    %d nonconvex e, f constraints added\n"%count)