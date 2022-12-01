import sys
import math
#from log import Logger
from gurobipy import *
import numpy as np
from myutils import breakexit
import time

epsilon4 = 1e-4

def lpformulator_ac(alldata):
  log = alldata['log']
  log.joint("AC formulation\n")

  starttime = time.time()

  retcode = lpformulator_setup(alldata)

  retcode, model, varcount = lpformulator_ac_body(alldata)


  model.params.NonConvex = 2

  model.params.DualReductions = 0

  #model.setParam(GRB.Param.MIPGap, 1.0e-10)
  #model.setParam(GRB.Param.FeasibilityTol, 1.0e-8)
  model.Params.MIPGap = 1.0e-3
  model.Params.OptimalityTol = 1.0e-3
  model.Params.FeasibilityTol = 1e-6 #1.0e-8

  feastol = model.Params.FeasibilityTol
  opttol = model.Params.OptimalityTol
  mipgap = model.Params.MIPGap
  log.joint('\n>> feastol %g opttol %g mipgap %g\n' %(feastol,opttol,mipgap))

  log.joint("variables = " + str(model.NumVars))
  log.joint("constraints = " + str(model.NumConstrs))

  model.optimize()

  if model.status == GRB.status.INF_OR_UNBD:
    log.joint('->LP infeasible or unbounded')
    model.Params.DualReductions = 0
    model.optimize()
  if model.status == GRB.status.INFEASIBLE:
    log.joint('->LP infeasible')                                            
    model.computeIIS()
    model.write("model.ilp")

  elif model.status == GRB.status.UNBOUNDED:
    log.joint('->LP unbounded')                                             

  elif model.status == GRB.OPTIMAL:
    log.joint(' ->Gurobi status OPTIMAL\n')


  log.joint('Optimal objective = %g' % model.objVal)

  model.printQuality()

  '''
  buses = alldata['buses']
  branches = alldata['branches']
  gens = alldata['gens']
  IDtoCountmap = alldata['IDtoCountmap']
  '''

def lpformulator_setup(alldata):
  retcode = 0

  alldata['maxdispersion_rad'] =   (math.pi/180.)*alldata['maxdispersion_deg']
  return retcode

def lpformulator_ac_body(alldata):
  retcode = varcount = 0

  log = alldata['log']
  model = Model('grbacs')
  
  retcode, varcount = lpformulator_ac_vars(alldata, model)

  if retcode == 0:
    log.joint('{} variables\n'.format(varcount))
    retcode = lpformulator_ac_constraints(alldata, model, varcount)

  log.joint('lpformulator_ac_body returns {}\n'.format(retcode))

  return retcode, model, varcount

def lpformulator_ac_vars(alldata, model):
  retcode = varcount = 0

  log = alldata['log']

  log.joint('creating variables\n')

  fixtolerance = 1e-05
  if alldata['fixtolerance'] > 0:
    fixtolerance = alldata['fixtolerance']

  numbuses = alldata['numbuses']
  buses = alldata['buses']
  IDtoCountmap = alldata['IDtoCountmap']
  
  #first, bus related variables

  cvar = {}
  svar = {}
  Pinjvar = {}
  Qinjvar = {}
  GenPvar = {}
  GenQvar = {}

  gens = alldata['gens']
  
  for j in range(1,numbuses+1):
    bus = buses[j]

    #first, injection variables
    maxprod = bus.Vmax*bus.Vmax
    minprod = bus.Vmin*bus.Vmin
          
    ubound = maxprod
    lbound = minprod

    if alldata['FIXCS'] and bus.inputvoltage:
      lbound = bus.inputV*bus.inputV - fixtolerance
      ubound = bus.inputV*bus.inputV + fixtolerance

    cvar[bus] = model.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "c_"+str(bus.nodeID)+"_"+str(bus.nodeID))

    bus.cffvarind = varcount
    varcount += 1

    # csdefslacks to be done

    Plbound = Qlbound = -GRB.INFINITY
    Pubound = Qubound = GRB.INFINITY

    Pubound, Plbound, Qubound, Qlbound = computebalbounds(log, alldata, bus)

    Pinjvar[bus] = model.addVar(obj = 0.0, lb = Plbound, ub = Pubound, name = "IP_"+str(bus.nodeID))
    bus.Pinjvarind = varcount
    varcount += 1
    Qinjvar[bus] = model.addVar(obj = 0.0, lb = Qlbound, ub = Qubound, name = "IQ_"+str(bus.nodeID))
    bus.Qinjvarind = varcount    
    varcount += 1

    #next, generator variables
    for genid in bus.genidsbycount:
      gen = gens[genid]
      lower = gen.Pmin*gen.status
      upper = gen.Pmax*gen.status
      #if bus.nodetype == 3:
      #  upper = GRB.INFINITY
      #  lower = -GRB.INFINITY  #ignoring slack bus
      GenPvar[gen] = model.addVar(obj = 0.0, lb = lower, ub = upper, name = "GP_"+str(gen.count)+"_"+str(gen.nodeID))
      gen.Pvarind = varcount
      varcount += 1
      lower = gen.Qmin*gen.status
      upper = gen.Qmax*gen.status
      if bus.nodetype == 3:
        upper = GRB.INFINITY
        lower = -GRB.INFINITY
      GenQvar[gen] = model.addVar(obj = 0.0, lb = lower, ub = upper, name = "GQ_"+str(gen.count)+"_"+str(gen.nodeID))
      gen.Qvarind = varcount
      varcount += 1

  alldata['LP']['GenPvar'] = GenPvar
  alldata['LP']['GenQvar'] = GenQvar  
  #next, branch related variables
  branches = alldata['branches']
  numbranches = alldata['numbranches']
  for j in range(1,1+numbranches):
    branch = branches[j]
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]
    busf = buses[count_of_f]
    bust = buses[count_of_t]
    maxprod = buses[count_of_f].Vmax*buses[count_of_t].Vmax
    minprod = buses[count_of_f].Vmin*buses[count_of_t].Vmin

    #Assumption 1.  zero angle difference is always allowed! More precisely minangle_rad <= 0 and maxaxangle_rad >= 0

    if branch.maxangle_rad < 0 or branch.minangle_rad > 0:
      log.joint('broken assumption 1: branch j {} f {} t {} minanglerad {} maxanglerad {}\n'.format(j, f, t, branch.minangle_rad, branch.maxangle_rad));

    
    ubound = ubasic = maxprod
    lbound = lbasic = -maxprod
    maxanglerad = branch.maxangle_rad
    minanglerad = branch.minangle_rad    

    #first, cosine

    if maxanglerad <=  .5*math.pi:
      #in this case minangle <= 0 
      if minanglerad >= -.5*math.pi:
        lbound = minprod*min(math.cos(maxanglerad), math.cos(minanglerad))
      elif minanglerad >= - math.pi:
        lbound = maxprod*math.cos(minangle_rad)  #which is negative
      elif minanglerad >= - 1.5*math.pi:
        lbound = -maxprod
      else:
        lbound = -maxprod
    elif maxanglerad <= math.pi:
      if minanglerad >= -.5*math.pi:
        lbound = maxprod*math.cos(maxanglerad)
      elif minanglerad >= - math.pi:
        lbound = maxprod*min(math.cos(maxanglerad), math.cos(minanglerad))
      elif minanglerad >= - 1.5*math.pi:
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
    
    cvar[branch] = model.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "c_"+str(j)+"_"+str(busf.nodeID) + "_" + str(bust.nodeID))

    branch.cftvarind = varcount
    varcount += 1

    #next, sine

    if maxanglerad <= math.pi/2:
      ubound = maxprod*sin(maxanglerad)

      if  minanglerad >= -.5*math.pi:
        lbound = maxprod*sin(minanglerad)
      elif  minanglerad >= - math.pi:
        lbound = -maxprod
      elif  minanglerad >= - 1.5*math.pi:
        ubound = maxprod*max( sin(maxanglerad), sin(minanglerad))
        lbound = -maxprod
      else:
        ubound = maxprod
        lbound = -maxprod 
    elif maxanglerad <= math.pi:
      ubound = maxprod

      if minanglerad >= -.5*math.pi:
        lbound = maxprod*sin(minanglerad)
      elif minanglerad >= - math.pi:
        lbound = -maxprod
      elif minanglerad >= - 1.5*math.pi:
        lbound = -maxprod
      else:
        lbound = -maxprod 

    elif maxanglerad <= 1.5*math.pi:
      ubound = maxprod

      if minanglerad >= -.5*math.pi:
        lbound = maxprod*min(sin(maxanglerad), sin(minanglerad))
      elif minanglerad >= - math.pi:
        lbound = -maxprod
      elif minanglerad >= - 1.5*math.pi:
        lbound = -maxprod
      else:
        lbound = -maxprod 
    else:
      ubound = maxprod
      lbound = -maxprod

    svar[branch] = model.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "s_"+str(j)+"_"+str(busf.nodeID) + "_" + str(bust.nodeID))

    branch.sftvarind = varcount
    varcount += 1

  #powerflow variables

  if alldata['use_ef'] and alldata['useconvexformulation'] == False:
    retcode,efvarcount = lpformulator_ac_add_efvars(alldata, model, varcount)
    if retcode:
      sys.exit("failure to add e, f, variables")
    varcount += efvarcount

  alldata['LP']['Pinjvar'] = Pinjvar
  alldata['LP']['Qinjvar'] = Qinjvar  

  Pvar_f = {}
  Qvar_f = {}
  Pvar_t = {}
  Qvar_t = {}

  for j in range(1,1+numbranches):
    branch = branches[j]
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]
    busf = buses[count_of_f]
    bust = buses[count_of_t]

    ubound = branch.limit
    lbound = -ubound

    Pvar_f[branch] = model.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "P_"+str(j)+"_"+str(busf.nodeID) + "_" + str(bust.nodeID))
    branch.Pftvarind = varcount
    varcount += 1

    Pvar_t[branch] = model.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "P_"+str(j)+"_"+str(bust.nodeID) + "_" + str(busf.nodeID))
    branch.Ptfvarind = varcount
    varcount += 1

    Qvar_f[branch] = model.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "Q_"+str(j)+"_"+str(busf.nodeID) + "_" + str(bust.nodeID))
    branch.Qftvarind = varcount
    varcount += 1
    
    Qvar_t[branch] = model.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "Q_"+str(j)+"_"+str(bust.nodeID) + "_" + str(busf.nodeID))
    branch.Qtfvarind = varcount
    varcount += 1
    
    
  lincostvar = model.addVar(obj = 1.0, lb = -GRB.INFINITY, ub = GRB.INFINITY, name = "lincost")
  alldata['LP']['lincostvar'] = lincostvar
  alldata['LP']['lincostvarind'] = varcount
  varcount += 1


  if alldata['usequadcostvar']:
    quadcostvar = model.addVar(obj = 1.0, lb = 0, ub = GRB.INFINITY, name = "quadcost")
    alldata['LP']['quadcostvar'] = quadcostvar
    alldata['LP']['quadcostvarind'] = varcount
    varcount += 1

  constobjval = 0  
  for gen in gens.values():  
    if gen.status > 0:
      constobjval += gen.costvector[ gen.costdegree ]
  
  constvar = model.addVar(obj = constobjval, lb = 1.0, ub = 1.0, name = "constant")
  alldata['LP']['constvar'] = constvar
  varcount += 1

  model.update()

  alldata['LP']['cvar'] = cvar
  alldata['LP']['svar'] = svar
  alldata['LP']['Pvar_f'] = Pvar_f
  alldata['LP']['Pvar_t'] = Pvar_t  
  alldata['LP']['Qvar_f'] = Qvar_f
  alldata['LP']['Qvar_t'] = Qvar_t  
  


  log.joint('lpformulator_ac_vars returns {}\n'.format(retcode))

  return retcode, varcount

def lpformulator_ac_add_efvars(alldata, model, varcount):
  retcode = efvarcount = 0
  
  log = alldata['log']

  log.joint('creating e,f variables\n')

  fixtolerance = 1e-05
  if alldata['fixtolerance'] > 0:
    fixtolerance = alldata['fixtolerance']

  numbuses = alldata['numbuses']
  buses = alldata['buses']
  IDtoCountmap = alldata['IDtoCountmap']

  evar = {}
  fvar = {}

  for j in range(1,numbuses+1):
    bus = buses[j]
    ubound = bus.Vmax
    lbound = -ubound

    if alldata['usemaxdispersion']:
      lbound = bus.Vmin*math.cos(alldata['maxdispersion_rad'])
      ubound = bus.Vmax

    evar[bus] = model.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "e_"+str(bus.nodeID))
    bus.evarind = varcount + efvarcount
    efvarcount += 1

    if alldata['usemaxdispersion']:
      lbound = 0
      ubound = bus.Vmax*math.sin(alldata['maxdispersion_rad'])

    elif j == alldata['refbus']:
      ubound = lbound = 0

      
    
    fvar[bus] = model.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "f_"+str(bus.nodeID))
    bus.fvarind = varcount + efvarcount
    efvarcount += 1
    
  alldata['LP']['evar'] = evar
  alldata['LP']['fvar'] = fvar  

  return retcode, efvarcount

def computebalbounds(log, all_data, bus):

  #first let's get max/min generations

  loud = 0

  baseMVA = all_data['baseMVA']
  gens = all_data['gens']

  Pubound = Plbound = 0
  Qubound = Qlbound = 0


  for gencounter in bus.genidsbycount:
    if gens[gencounter].status:
      Pubound += gens[gencounter].Pmax
      Plbound += gens[gencounter].Pmin
      Qubound += gens[gencounter].Qmax
      Qlbound += gens[gencounter].Qmin

    if loud:
     #log.joint(" Pubound for " + str(bus.nodeID) + " " + str(Pubound) + " genc " + str(gencounter) + "\n")
     #log.joint(" Plbound for " + str(bus.nodeID) + " " + str(Plbound) + " genc " + str(gencounter) + "\n")
     log.joint(" Qubound for " + str(bus.nodeID) + " " + str(Qubound) + " genc " + str(gencounter) + "\n")
     log.joint(" Qlbound for " + str(bus.nodeID) + " " + str(Qlbound) + " genc " + str(gencounter) + "\n")

  Pubound -= bus.Pd
  Plbound -= bus.Pd
  Qubound -= bus.Qd
  Qlbound -= bus.Qd

  if bus.nodetype == 4:
    Pubound = Plbound = Qubound = Qlbound = 0

  if loud:
     #log.joint(" Pubound for " + str(bus.nodeID) + " final " + str(Pubound) + "\n")
     #log.joint(" (Pd was %g)\n" %bus.Pd)
     #log.joint(" Plbound for " + str(bus.nodeID) + " final " + str(Plbound) + "\n")
     log.joint(" Qubound for " + str(bus.nodeID) + " final " + str(Qubound) + "\n")
     log.joint(" (Qd was %g)\n" %bus.Qd)
     log.joint(" Qlbound for " + str(bus.nodeID) + " final " + str(Qlbound) + "\n")
     breakexit(" ")
  
  return Pubound, Plbound, Qubound, Qlbound

def lpformulator_ac_constraints(alldata, model, varcount):
  retcode = 0

  log = alldata['log']
  log.joint('building constraints\n')
    
  numbuses    = alldata['numbuses']
  buses = alldata['buses']
  numbranches = alldata['numbranches']
  branches = alldata['branches']
  gens = alldata['gens']
  IDtoCountmap = alldata['IDtoCountmap']
  cvar = alldata['LP']['cvar']
  svar = alldata['LP']['svar']
  Pvar_f = alldata['LP']['Pvar_f']
  Pvar_t = alldata['LP']['Pvar_t']    
  Qvar_f = alldata['LP']['Qvar_f']
  Qvar_t = alldata['LP']['Qvar_t']
  Pinjvar = alldata['LP']['Pinjvar']
  Qinjvar = alldata['LP']['Qinjvar']  


  print("cost def constraints\n")
  
  GenPvar = alldata['LP']['GenPvar']
  GenQvar = alldata['LP']['GenQvar']  
  lincostvar = alldata['LP']['lincostvar']


  model.update()

  coeff = [gen.costvector[ gen.costdegree - 1] for gen in gens.values()]
  variables = [GenPvar[gen] for gen in gens.values()]
  expr = LinExpr(coeff, variables)
  model.addConstr(expr == lincostvar, name="lincostdef")

  numquadgens = 0
  for gen in gens.values():
    '''
    print(gen.count, gen.nodeID)
    print(gen.costdegree, gen.costvector)
    '''
    if gen.costdegree >= 2 and gen.costvector[0] > 0 and gen.status:
      numquadgens += 1
  log.joint("number of generators with quadratic cost coefficient: {}\n".format(numquadgens))

  if numquadgens > 0:
    if alldata['usequadcostvar']:
      quadcostvar = alldata['LP']['quadcostvar']
      log.joint('adding quadcost def constraint\n')
      qcost = QuadExpr()
      for gen in gens.values():
        if gen.costdegree == 2 and gen.costvector[0] != 0:
          qcost += gen.costvector[0]*GenPvar[gen]*GenPvar[gen]
      model.addConstr(qcost <= quadcostvar, name = "qcostdef")
    else:
      log.joint('adding quad cost to objective\n')
      oldobj = model.getObjective()
      for gen in gens.values():
        if gen.costdegree == 2 and gen.costvector[0] != 0:
          oldobj += gen.costvector[0]*GenPvar[gen]*GenPvar[gen]
      model.setObjective(oldobj, GRB.MINIMIZE)
      

  #define flow variables
  log.joint('active power flow defs\n')
  for j in range(1,1+numbranches):
    branch = branches[j]
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]
    busf = buses[count_of_f]
    bust = buses[count_of_t]

    if alldata['substitute_nonconv']==False or alldata['use_ef'] == False:

      #  Gff cff + Gft cft + Bft sft
      constrname = "Pdef_"+str(j)+"_"+str(f)+"_"+str(t)
      expr = LinExpr()
      expr += branch.Gff*cvar[busf]
      expr += branch.Gft*cvar[branch]
      expr += branch.Bft*svar[branch]
      model.addConstr(expr == Pvar_f[branch], name = constrname)

      #  Gtt ctt + Gtf cft + Btf stf = Gtt ctt + Gtf cft - Btf sft
      constrname = "Pdef_"+str(j)+"_"+str(t)+"_"+str(f)
      expr = LinExpr()
      expr += branch.Gtt*cvar[bust]
      expr += branch.Gtf*cvar[branch]
      expr += -branch.Btf*svar[branch] #minus because svarft = -svartf
      model.addConstr(expr == Pvar_t[branch], name = constrname)
    else:
      breakexit("se")
  

  log.joint('reactive power flow defs\n')
  for j in range(1,1+numbranches):
    branch = branches[j]
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]
    busf = buses[count_of_f]
    bust = buses[count_of_t]

    if alldata['substitute_nonconv']==False or alldata['use_ef'] == False:
      # -Bff cff - Bft cft + Gft sft
      constrname = "Qdef_"+str(j)+"_"+str(f)+"_"+str(t)
      expr = LinExpr()
      expr += -branch.Bff*cvar[busf]
      expr += -branch.Bft*cvar[branch]
      expr += +branch.Gft*svar[branch]
      model.addConstr(expr == Qvar_f[branch], name = constrname)

      #  -Btt ctt - Btf cft + Gtf stf = -Btt ctt - Btf cft - Gtf sft 
      constrname = "Qdef_"+str(j)+"_"+str(t)+"_"+str(f)
      expr = LinExpr()
      expr += -branch.Btt*cvar[bust]
      expr += -branch.Btf*cvar[branch]
      expr += -branch.Gtf*svar[branch] #again, same minus
      model.addConstr(expr == Qvar_t[branch], name = constrname)
      
    else:
      breakexit('se')

  log.joint('balance constraints\n')

  for j in range(1,1+numbuses):
    bus = buses[j]
    constrname = "PBaldef"+str(j)+'_'+str(bus.nodeID)
    expr = LinExpr()
    for branchid in bus.frombranchids.values():
      expr += Pvar_f[ branches[branchid] ]
    for branchid in bus.tobranchids.values():
      expr += Pvar_t[ branches[branchid] ]


    if bus.Gs != 0:
      expr += bus.Gs*cvar[bus]
      
      
    model.addConstr(expr == Pinjvar[bus], name = constrname)

  for j in range(1,1+numbuses):
    bus = buses[j]
    constrname = "QBaldef"+str(j)+'_'+str(bus.nodeID)
    expr = LinExpr()
    for branchid in bus.frombranchids.values():
      expr += Qvar_f[ branches[branchid] ]
    for branchid in bus.tobranchids.values():
      expr += Qvar_t[ branches[branchid] ]
 
    if bus.Bs != 0:
      expr += (-bus.Bs)*cvar[bus]

    model.addConstr(expr == Qinjvar[bus], name = constrname)


  log.joint('injection definition constraints\n')
  
  for j in range(1,1+numbuses):
    bus = buses[j]
    constrname = "Bus_PInj_"+str(j)
    expr = LinExpr()

    if len(bus.genidsbycount) > 0:
      for genid in bus.genidsbycount:
        gen = gens[genid]
        expr += GenPvar[gen]

    model.addConstr(Pinjvar[bus] == expr - bus.Pd, name = constrname)

    constrname = "Bus_QInj_"+str(j)
    expr = LinExpr()

    if len(bus.genidsbycount) > 0:
      for genid in bus.genidsbycount:
        gen = gens[genid]
        expr += GenQvar[gen]

    model.addConstr(Qinjvar[bus] == expr - bus.Qd, name = constrname)

  log.joint('branch limits\n')
  for j in range(1,1+numbranches):
    branch = branches[j]
    if branch.status and branch.unboundedlimit == False:
      f = branch.f
      t = branch.t
      count_of_f = IDtoCountmap[f]
      count_of_t = IDtoCountmap[t]
      busf = buses[count_of_f]
      bust = buses[count_of_t]
      constrname = "limit_f_"+str(j)+"_"+str(f)+"_"+str(t)
      limexp = QuadExpr()
      limexp += Pvar_f[branch]*Pvar_f[branch] + Qvar_f[branch]*Qvar_f[branch]
      model.addConstr(limexp <= branch.limit**2, name = constrname)

      constrname = "limit_t_"+str(j)+"_"+str(t)+"_"+str(f)
      limexp = QuadExpr()
      limexp += Pvar_t[branch]*Pvar_t[branch] + Qvar_t[branch]*Qvar_t[branch]
      #themodel.cbLazy(limexp <= branch.limit**2)
      model.addConstr(limexp <= branch.limit**2, name = constrname)


  if alldata['skipjabr'] == False:
    log.joint('Jabr constraints\n')
    for j in range(1,1+numbranches):
      branch = branches[j]
      if branch.status:
        f = branch.f
        t = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]
        busf = buses[count_of_f]
        bust = buses[count_of_t]
        constrname = 'jabr_'+str(j)+'_'+str(f)+'_'+str(t)
        limexp = QuadExpr()
        limexp += cvar[branch]*cvar[branch] + svar[branch]*svar[branch]
        model.addConstr(limexp <= cvar[busf]*cvar[bust], name = constrname)
      else:
        log.joint('skipping Jabr inequalities\n')  

  if alldata['use_ef'] and alldata['useconvexformulation']==False:
    lpformulator_ac_add_nonconvexconstraints(alldata, model)
    
  model.update()
  

  model.write('foo.lp')  

  log.joint('lpformulator_ac_constraints returns {}\n'.format(retcode))

  return retcode

def lpformulator_ac_add_nonconvexconstraints(alldata,model):
  log = alldata['log']
  log.joint('adding nonconvex constraints\n')

  buses = alldata['buses']
  numbuses = alldata['numbuses']
  branches = alldata['branches']
  numbranches = alldata['numbranches']
  evar = alldata['LP']['evar']
  fvar = alldata['LP']['fvar']
  cvar = alldata['LP']['cvar']
  svar = alldata['LP']['svar']    
  IDtoCountmap = alldata['IDtoCountmap']
  
  log.joint('e,f nonconvex constraints\n')
  for j in range(1,1+numbuses):
    bus = buses[j]
    constrname = 'cbusdef_'+str(j)+'_'+str(bus.nodeID)
    model.addConstr(-cvar[bus] + evar[bus]*evar[bus] + fvar[bus]*fvar[bus] == 0, name =constrname)



  for j in range(1,1+numbranches):
    branch = branches[j]
    if branch.status:
      f = branch.f
      t = branch.t
      count_of_f = IDtoCountmap[f]
      count_of_t = IDtoCountmap[t]
      
      busf = buses[count_of_f]
      bust = buses[count_of_t]

      constrname = 'cdef_'+str(j)+'_'+str(f)+'_'+str(t)
      model.addConstr(-cvar[branch] + evar[busf]*evar[bust] + fvar[busf]*fvar[bust] == 0, constrname)
      constrname = 'sdef_'+str(j)+'_'+str(f)+'_'+str(t)
      model.addConstr(-svar[branch] - evar[busf]*fvar[bust] + fvar[busf]*evar[bust] == 0, constrname)
