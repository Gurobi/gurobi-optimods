import math
import time

import gurobipy as gp
from gurobipy import GRB
from myutils import break_exit
from grbfile import grbreadvoltsfile
from grbgraphical import grbgraphical
from grbformulator_ac import computebalbounds
import numpy as np

globalalldata = {}



def lpformulator_dc(alldata):
    """Formulate DCOPF model and solve it"""

    log = alldata['log']
    log.joint("\nDC formulation.\n")

    starttime = time.time()

    global globalalldata 
    globalalldata = alldata


    # Handle special settings
    #lpformulator_setup(alldata)

    sol_count = 0
    # Create model
    with gp.Env() as env, gp.Model('grbacs', env=env) as model:
        # Add model variables and constraints
        lpformulator_dc_body(alldata, model)

        if alldata['strictcheckvoltagesolution']:
            #check input solution against formulation
            spitoutvector = True
            #feascode = lpformulator_dc_strictchecker(alldata, model, spitoutvector) #uncomment this for check

        sol_count = lpformulator_dc_opt(alldata, model)

        endtime = time.time()
        log.joint("Overall time taken (model construction + optimization): %f s.\n"%(endtime-starttime))
        log.joint('Solution count: %d\n'%(sol_count))

        break_exit('Solver returned.')

        if sol_count > 0:
            lpformulator_dc_examine_solution(alldata, model)
    '''
    buses        = alldata['buses']
    branches     = alldata['branches']
    gens         = alldata['gens']
    IDtoCountmap = alldata['IDtoCountmap']
    '''

def lpformulator_dc_body(alldata, model):
    """Helper function for adding variables and constraints to the model"""

    log = alldata['log']

    # Create model variables
    lpformulator_dc_create_vars(alldata, model)
    # Create model constraints
    lpformulator_dc_create_constraints(alldata, model)


    model.update() # Update to get correct model stats
    log.joint("Constructed DCOPF model with %d variables and %d constraints.\n\n"%(model.NumVars, model.NumConstrs))

    model.write(alldata['lpfilename'])  # FIXME remove.  Jarek: I am using this for debugging, for now
    log.joint('Wrote LP to ' + alldata['lpfilename'] + '\n')
    #break_exit('wrote lp') # 

    alldata['model'] = model

def lpformulator_dc_create_vars(alldata, model):
    """Create model variables for DCOPF"""

    log = alldata['log']
    log.joint('Creating variables.\n')

    fixtolerance = 1e-05
    if alldata['fixtolerance'] > 0:
        fixtolerance = alldata['fixtolerance']

    numbuses     = alldata['numbuses']
    buses        = alldata['buses']
    IDtoCountmap = alldata['IDtoCountmap']
    gens     = alldata['gens']
    
    varcount = 0
    thetavar     = {}
    Pinjvar      = {}
    Pvar_f       = {}   #DC, so f-flow = t-flow
    twinPvar_f   = {}   #auxiliary variable in case branch-switching is being used
    GenPvar      = {}

    for j in range(1,numbuses+1):
        bus       = buses[j]
        ubound = 2*math.pi
        lbound = -ubound

        if bus.inputvoltage:
            candidatelbound = bus.inputA_rad - fixtolerance
            candidateubound = bus.inputA_rad + fixtolerance
            lbound = max(lbound, candidatelbound)
            ubound = min(ubound, candidateubound)
        
        thetavar[bus]  = model.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "theta_"+str(bus.nodeID))
        bus.thetavarind = varcount
        varcount += 1
        
        Plbound = Qlbound = -GRB.INFINITY
        Pubound = Qubound = GRB.INFINITY
        Pubound, Plbound, Qubound, Qlbound = computebalbounds(log, alldata, bus)

        Pinjvar[bus] = model.addVar(obj = 0.0, lb = Plbound, ub = Pubound,
                                    name = "IP_%d"%bus.nodeID)
        bus.Pinjvarind = varcount
        #comment: Pinjvar is the variable modeling total active power injected by bus j into the branches incident with j
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

        
    

    alldata['LP']['thetavar'] = thetavar
    alldata['LP']['Pinjvar']  = Pinjvar
    alldata['LP']['GenPvar']  = GenPvar    

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
        if branch.constrainedflow:
            ubound     = branch.limit
        else:
            ubound     = alldata['sumPd'] #DC
        lbound     = -ubound

        Pvar_f[branch] = model.addVar(obj = 0.0, lb = lbound, ub = ubound,
                                      name = "P_%d_%d_%d"%(j, busf.nodeID, bust.nodeID))
        branch.Pftvarind = varcount
        varcount += 1

        if alldata['branchswitching_mip']:
            twinPvar_f[branch] = model.addVar(obj = 0.0, lb = lbound, ub = ubound,
                                      name = "twinP_%d_%d_%d"%(j, busf.nodeID, bust.nodeID))
            branch.twinPftvarind = varcount
            varcount += 1


    alldata['LP']['Pvar_f']     = Pvar_f
    alldata['LP']['twinPvar_f'] = twinPvar_f    

    zvar = {}
    if alldata['branchswitching_mip']:
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
    
def lpformulator_dc_create_constraints(alldata, model):
    """"Create constraint for ACOPF"""

    numbuses     = alldata['numbuses']
    buses        = alldata['buses']
    numbranches  = alldata['numbranches']
    branches     = alldata['branches']
    gens         = alldata['gens']
    IDtoCountmap = alldata['IDtoCountmap']
    log          = alldata['log']

    thetavar     = alldata['LP']['thetavar']
    Pvar_f       = alldata['LP']['Pvar_f']
    twinPvar_f   = alldata['LP']['twinPvar_f']
    zvar         = alldata['LP']['zvar']
    Pinjvar      = alldata['LP']['Pinjvar']
    GenPvar      = alldata['LP']['GenPvar']
    lincostvar   = alldata['LP']['lincostvar']
    log          = alldata['log']


    log.joint("Creating constraints.\n")
    log.joint("  Adding cost definition.\n")

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
            log.joint("    Adding quad cost to objective.\n")
            model.update() # Necessary to flush changes in the objective function
            oldobj = model.getObjective()
            newobj = gp.QuadExpr(oldobj)
            for gen in gens.values():
                if gen.costdegree == 2 and gen.costvector[0] != 0:
                    newobj.add(gen.costvector[0]*GenPvar[gen]*GenPvar[gen])

            model.setObjective(newobj, GRB.MINIMIZE)

    # Active PF defs
    log.joint("  Adding active power flow definitions.\n")
    count = 0
    for j in range(1,1+numbranches):
        branch     = branches[j]
        f          = branch.f
        t          = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]
        busf       = buses[count_of_f]
        bust       = buses[count_of_t]
        branch.Pfcname = "Pdef_%d_%d_%d"%(j, f, t)
        if branch.status:
            # Pf = (thetaf - thetat)/(x*ratio)
            coeff = 1/(branch.x*branch.ratio)
            exp = Pvar_f[branch]
            if alldata['branchswitching_mip']:
                exp += twinPvar_f[branch]
            #angle_exp = coeff*thetavar[busf] - coeff*thetavar[bust] - coeff*branch.angle_rad
            branch.Pdeffconstr = model.addConstr(exp == coeff*thetavar[busf] - coeff*thetavar[bust] - coeff*branch.angle_rad, name = branch.Pfcname)
            count += 1
            #print(count, coeff, branch.angle_rad, angle_exp)
            #break_exit('pdff')

            if alldata['branchswitching_mip']:
                coeff = branch.limit
                model.addConstr(Pvar_f[branch] <=  coeff*zvar[branch], name = 'upmip_%d_%d_%d'%(j, f, t))
                model.addConstr(Pvar_f[branch] >= -coeff*zvar[branch], name = 'dnmip_%d_%d_%d'%(j, f, t))
                model.addConstr(twinPvar_f[branch] <=  coeff*(1-zvar[branch]), name = 'upmip_twin_%d_%d_%d'%(j, f, t))
                model.addConstr(twinPvar_f[branch] >= -coeff*(1 - zvar[branch]), name = 'dnmip_twin_%d_%d_%d'%(j, f, t))
        else:
            branch.Pdeffconstr = model.addConstr( Pvar_f[branch] == 0, name = branch.Pfcname)


    log.joint("    %d active power flow definitions added.\n"%count)
                  

    # Balance constraints
    log.joint("  Adding constraints stating bus injection = total outgoing power flow.\n")
    count = 0
    for j in range(1,1+numbuses):
        bus  = buses[j]
        expr = gp.LinExpr()
        for branchid in bus.frombranchids.values():
            expr.add(Pvar_f[branches[branchid]])

        for branchid in bus.tobranchids.values():
            expr.add(-Pvar_f[branches[branchid]])

        model.addConstr(expr == Pinjvar[bus], name = "PBaldef%d_%d"%(j, bus.nodeID))            
        
        count += 1
    log.joint("    %d constraints added.\n"%count)

    # Injection defs
    log.joint("  Adding injection definition constraints.\n")
    count = 0
    for j in range(1,1+numbuses):
        bus  = buses[j]
        expr = gp.LinExpr()

        if len(bus.genidsbycount) > 0:
            for genid in bus.genidsbycount:
                gen = gens[genid]
                expr.add(GenPvar[gen])

        model.addConstr(Pinjvar[bus] == expr - bus.Pd, name = "Bus_PInj_%d"%j)
        count += 1

    log.joint("    %d injection definition constraints added.n"%count)

    if alldata['branchswitching_mip']:
        boundzs = True
        if boundzs:
            exp = gp.LinExpr()
            delta = 50
            N = numbranches - delta  # <<<<<<---- here is the heuristic lower bound
            for j in range(1,1+numbranches):
                branch     = branches[j]
                exp += zvar[branch]
            model.addConstr(exp >= N, name = "sumzbd")


def lpformulator_dc_opt(alldata, model):

    log = alldata['log']
    numbranches  = alldata['numbranches']
    branches     = alldata['branches']

    
    # Specific settings for better convergence
    #model.params.DualReductions = 0
    
    #model.setParam(GRB.Param.MIPGap, 1.0e-10)
    #model.setParam(GRB.Param.FeasibilityTol, 1.0e-8)
    model.Params.MIPGap         = 1.0e-2
    model.Params.OptimalityTol  = 1.0e-2
    model.Params.FeasibilityTol = 1.0e-6 # 1.0e-8
    
    feastol = model.Params.FeasibilityTol
    opttol  = model.Params.OptimalityTol
    mipgap  = model.Params.MIPGap

    model.Params.SolutionLimit = 10

    if alldata['branchswitching_mip']:
        log.joint('Using mip start with all branches kept on\n')
        #mip start
        zvar = alldata['LP']['zvar']
        for j in range(1,1+numbranches):
            branch = branches[j]
            zvar[branch].Start = 1.0
        

        zholder = np.zeros(numbranches)
        alldata['LP']['zholder'] = zholder
        
        # Optimize
        model._vars = model.getVars()
        model.optimize(mycallback)
    else:
        model.optimize()
        
    # Check model status and re-optimize or try computing an IIS if necessary
    if model.status == GRB.INF_OR_UNBD:
        log.joint("\nModel Status: infeasible or unbounded.\n")
        log.joint("Re-optimizing with DualReductions turned off.\n\n")
        model.Params.DualReductions = 0
        model.optimize()
        
    if model.status == GRB.INFEASIBLE:
        log.joint("\nModel Status: infeasible.")
        log.joint("Computing IIS...\n\n")
        model.computeIIS()
        log.joint("\nIIS computed, writing IIS to file acopfmodel.ilp\n\n")
        model.write("acopfmodel.ilp")
        
    elif model.status == GRB.UNBOUNDED:
        log.joint("\nModel Status: unbounded.\n\n")

    elif model.status == GRB.INTERRUPTED:
        log.joint("\nModel Status: interrupted.\n\n")

    elif model.status == GRB.OPTIMAL:
        log.joint("\nModel Status: optimal.\n\n")
        
    # Only print objective value and solution quality if at least
    # one feasible point is available
    if model.SolCount > 0:
            log.joint("Objective value = %g"%model.objVal)
            model.printQuality()
            # FIXME yes, to be done
            # Here we should gather optimal solution values and gather them
            # in a standardized format which ultimately will be returned to the user
            # 

    return model.SolCount


def lpformulator_dc_examine_solution(alldata, model):
    log = alldata['log']

    numbuses     = alldata['numbuses']
    buses        = alldata['buses']
    branches     = alldata['branches']
    numbranches  = alldata['numbranches']
    branches     = alldata['branches']
    gens         = alldata['gens']
    IDtoCountmap = alldata['IDtoCountmap']
    log          = alldata['log']
    if alldata['branchswitching_mip']: zholder = alldata['LP']['zholder']

    thetavar     = alldata['LP']['thetavar']
    Pvar_f       = alldata['LP']['Pvar_f']
    twinPvar_f   = alldata['LP']['twinPvar_f']
    zvar         = alldata['LP']['zvar']
    Pinjvar      = alldata['LP']['Pinjvar']
    GenPvar      = alldata['LP']['GenPvar']
    lincostvar   = alldata['LP']['lincostvar']

    #for var in model.getVars():
    #    print(var.varname)

    loud = False
    numzeros = 0
    for j in range(1,1+numbranches):
        branch     = branches[j]
        f          = branch.f
        t          = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]
        busf       = buses[count_of_f]
        bust       = buses[count_of_t]

        #print(Pvar_f[branch].x)

        if alldata['branchswitching_mip']:
            zholder[j-1] = zvar[branch].x
            if zvar[branch].x < 0.5:  #turned off
                numzeros += 1
                if loud:
                    log.joint('branch %d (%d, %d) switched off\n'%(j, f, t))

        #log.joint('numzeros: %d\n'%numzeros)
                
    log.joint(' Done examining solution.\n')
    if alldata['dographics']:
        grbgraphical(alldata, 'branchswitching')
    
def mycallback(model, where):
    if where == GRB.Callback.MIPSOL:
        #solcnt = model.cbGet(GRB.Callback.MIP_SOLCNT)
        #print(solcnt)
        x = model.cbGetSolution(model._vars)
        log = globalalldata['log']
        log.joint('Inside callback.\n')
        objval = model.cbGet(GRB.Callback.MIPSOL_OBJ)
        log.joint('Found new solution of value %.3e\n'%objval)

        numbuses     = globalalldata['numbuses']
        buses        = globalalldata['buses']
        numbranches  = globalalldata['numbranches']
        branches     = globalalldata['branches']
        gens         = globalalldata['gens']
        IDtoCountmap = globalalldata['IDtoCountmap']
        log          = globalalldata['log']

        thetavar     = globalalldata['LP']['thetavar']
        Pvar_f       = globalalldata['LP']['Pvar_f']
        twinPvar_f   = globalalldata['LP']['twinPvar_f']
        zvar         = globalalldata['LP']['zvar']
        zholder = globalalldata['LP']['zholder']
        numzeros = 0        
        for j in range(1, 1+numbranches):
            branch = branches[j]
            zholder[j-1] = x[branch.switchvarind]
            if x[branch.switchvarind] < 0.5:
                numzeros += 1
                #print(j, x[branch.switchvarind])

        #log.joint('numzeros: %d\n'%numzeros)
        #break_exit('HEY HEY HEY')
        grbgraphical(globalalldata, 'branchswitching')
        #break_exit('HO HO HO')
