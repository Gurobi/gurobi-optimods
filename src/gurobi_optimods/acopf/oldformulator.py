import sys
import math
#from log import Logger
from gurobipy import *
import numpy as np
from myutils import breakexit
from digitizer import *
import reader
import time

epsilon4 = 1e-4
  
def lpformulator(log, all_data):
  log.joint("formulating in memory\n")
  if all_data['doac']:
    lpformulator_ac(log, all_data)
  elif all_data['dodcbasic']:
    lpformulator_dc_basic(log, all_data)
  elif all_data['dodc']:
    lpformulator_dc(log, all_data)


def lpformulator_dc_basic(log, all_data):
  log.joint("DC formulation\n")

  buses = all_data['buses']
  branches = all_data['branches']
  gens = all_data['gens']
  baseMVA = all_data['baseMVA']
  IDtoCountmap = all_data['IDtoCountmap']
  
  pi = math.pi

  themodel = Model("dcmodel")

  thetavar = {}
  Pvar = {}
  InjPvar = {}
  GenPvar = {}
  
  epsilon3 = all_data['epsilon_3']

  for bus in buses.values():

    myub = pi/2
    mylb = -pi/2


    if all_data['voltsfilename'] != 'NONE' and bus.nodeID in all_data['inputvolts']:
      myub = all_data['inputvolts'][bus.nodeID][1] + epsilon3
      mylb = all_data['inputvolts'][bus.nodeID][1] - epsilon3
   
    if bus.nodeID == all_data['slackbus']:
      mylb = myub = 0

    thetavar[bus] = themodel.addVar(obj = 0.0, lb = mylb, ub = myub, name = "theta_"+str(bus.nodeID))

  log.joint("angle variables\n")    

  sumlower = 0
  for bus in buses.values():
    Pubound, Plbound, Qubound, Qlbound = computebalbounds(log, all_data, bus)

    # net sum of flows on branches = generation   - demand
    # qtty on the right is what we represent using the variable InjPvar
    # Pubound and Plbound are the upper and lower bounds on generation - demand

    if bus.nodetype == 3:
      Pubound = GRB.INFINITY
      Plbound = -GRB.INFINITY
    else:
      sumlower += Plbound

    
    InjPvar[bus] = themodel.addVar(obj = 0.0, lb = Plbound, ub = Pubound, name = "IP_"+str(bus.nodeID))
    
  log.joint("balance variables\n")
  
  for bus in buses.values():

    for genid in bus.genidsbycount:
      gen = gens[genid]
      lower = gen.Pmin*gen.status
      upper = gen.Pmax*gen.status
  
      if bus.nodetype == 3:
        upper = GRB.INFINITY
        lower = -GRB.INFINITY
      GenPvar[gen] = themodel.addVar(obj = 0.0, lb = lower, ub = upper, name = "GP_"+str(gen.count)+"_"+str(gen.nodeID))
      
  log.joint("generator variables\n")

  for branch in branches.values():
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]
    #Pvar
    ubound = branch.limit
    lbound = -branch.limit
    if all_data['flowsfilename'] != 'NONE' and branch.count in all_data['inputPf']:
      ubound = all_data['inputPf'][branch.count] + epsilon3
      lbound = all_data['inputPf'][branch.count] - epsilon3

    Pvar[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "P_"+str(branch.count)+","+str(f) + "," + str(t))
  log.joint("branch variables\n")

  lincostvar = themodel.addVar(obj = 1.0, lb = -GRB.INFINITY, ub = GRB.INFINITY, name = "lincost")
  quadcostvar = themodel.addVar(obj = 1.0, lb = 0, ub = GRB.INFINITY, name = "quadcost")

  constobjval = 0
  for gen in gens.values():
    if gen.status > 0:
      constobjval += gen.costvector[2]

  # if the bus also participates in correction we adjust the constant objective alue

  log.joint(" constant in objective: " + str(constobjval) + "\n")
    
  
  constvar = themodel.addVar(obj = constobjval, lb = 1.0, ub = 1.0, name = "constant")

  themodel.update()

  #constraints

  #cost definitions
  
  qcost = QuadExpr()
  
  for gen in gens.values():
    if gen.costdegree==2 and gen.costvector[0] != 0:
      qcost += gen.costvector[0]*GenPvar[gen]*GenPvar[gen]
  
  themodel.addConstr(qcost <= quadcostvar, name = "qcostdef")

  coeff = [gen.costvector[ gen.costdegree - 1] for gen in gens.values()]
  variables = [GenPvar[gen] for gen in gens.values()]
  expr = LinExpr(coeff, variables)
  themodel.addConstr(expr == lincostvar, name="lincostdef")

  log.joint("cost def constraints\n")

  #define flow variables

  for branch in branches.values():
    f = branch.f
    t = branch.t
    coeff = 1/(branch.x*branch.ratio)
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]

    constrname = "Pdef_"+str(branch.count)+","+str(f)+","+str(t)
    expr = LinExpr()
    expr += coeff*thetavar[buses[count_of_f]]
    expr += -coeff*thetavar[buses[count_of_t]]
    themodel.addConstr(Pvar[branch] == expr - coeff*branch.angle_rad, name = constrname)

  log.joint("flow var def constraints\n")

  
  for bus in buses.values():
    nodeID = bus.nodeID
    constrname = "PBaldef"+str(nodeID)
    outflowexpr = LinExpr()
    for branchid in bus.frombranchids.values():
      outflowexpr += Pvar[ branches[branchid] ]
    for branchid in bus.tobranchids.values():
      outflowexpr -= Pvar[ branches[branchid] ]
      
    themodel.addConstr(outflowexpr == InjPvar[bus], name = constrname)

  log.joint("balance constraints\n")

  for bus in buses.values():

    constrname = "Bus_PInjdef_"+str(bus.nodeID)
    expr = LinExpr()

    if len(bus.genidsbycount) > 0:
      for genid in bus.genidsbycount:
        gen = gens[genid]
        expr += GenPvar[gen]

    themodel.addConstr(InjPvar[bus] == expr - bus.Pd, name = constrname)

  log.joint("generator def constraints\n")      

  themodel.update()

  log.joint("writing to lpfile " + all_data['lpfilename'] + "\n")
  themodel.write(all_data['lpfilename'])
  themodel.optimize()

  if themodel.status == GRB.status.INF_OR_UNBD:
    log.joint('->LP infeasible or unbounded\n')
    return themodel.status, 0

  if themodel.status == GRB.status.INFEASIBLE:
    log.joint('->LP infeasible\n')
    return themodel.status, 0

  if themodel.status == GRB.status.UNBOUNDED:
    log.joint('->LP unbounded\n')
    return themodel.status, 0

  log.joint(' objective = %g\n' % themodel.objVal)
  loud = 0
  totalqcost = totallcost = 0
  for bus in buses.values():
    if len(bus.genidsbycount) > 0:
      for genid in bus.genidsbycount:
        gen = gens[genid]
        genvalue = GenPvar[gen].x
        lincosterm = gen.costvector[1]*genvalue
        quadcosterm = gen.costvector[0]*genvalue*genvalue
        costerm = lincosterm + quadcosterm
        totalqcost += quadcosterm
        totallcost += lincosterm
        if loud:
          log.joint(" bus " + str(bus.nodeID) + " gen " + str(genid) + " -> " + str(genvalue) + " qcost " + str(quadcosterm) + " lcost " + str(lincosterm) + "\n")

  log.joint(" total qcost: " + str(totalqcost) + "\n")
  log.joint(" total lcost: " + str(totallcost) + "\n")
  totalvariablecost = totalqcost + totallcost
  log.joint(" total variable cost: " + str(totalvariablecost) + "\n")

  
  loud = 0

  psol = {}
  for geni in gens.values():
    generationi = GenPvar[geni].x
    #print geni.nodeID, generationi
    psol[geni] = generationi

  costestimator = 0
  for geni in gens.values():
    i = geni.count
    c0i = geni.costvector[0]
    c1i = geni.costvector[1]
    pi = psol[geni]

    log.joint("generator " + str(i) + " c0 " + str( c0i) + " c1 " + str(c1i) +" p " + str(pi) + " min " + str(geni.Pmin*gen.status) + " max " + str(geni.Pmax*gen.status) + "\n")

    costestimator += c0i*(pi*pi) + c1i*pi

  log.joint(" costestimator = " + str(costestimator) + "\n")

#  for branch in branches.values():
#    f = branch.f
#    t = branch.t
#    idbranch = branch.count
#    log.joint(str(branch.count) + " " + str(f) + " " + str(t) + " " + str(Pvar[branch].x*baseMVA) + "\n")
  
    

def lpformulator_dc(log, all_data):
  log.joint("DC formulation\n")

  buses = all_data['buses']
  branches = all_data['branches']
  gens = all_data['gens']
  IDtoCountmap = all_data['IDtoCountmap']
  stocs = all_data['stocs']
  alphas = all_data['alphas']

  
  pi = math.pi

  themodel = Model("dcmodel")

  thetavar = {}
  Pvar = {}
  InjPvar = {}
  GenPvar = {}
  
  epsilon3 = 1e-2

  for bus in buses.values():

    myub = len(buses)*pi
    mylb = -myub


    if all_data['voltsfilename'] != 'NONE' and bus.nodeID in all_data['inputvolts']:
      myub = all_data['inputvolts'][bus.nodeID][1] + epsilon3
      mylb = all_data['inputvolts'][bus.nodeID][1] - epsilon3
   
    if bus.nodeID == all_data['slackbus']:
      mylb = myub = 0

    thetavar[bus] = themodel.addVar(obj = 0.0, lb = mylb, ub = myub, name = "theta_"+str(bus.nodeID))

  log.joint("angle variables\n")    

  sumlower = 0
  for bus in buses.values():
    Pubound, Plbound, Qubound, Qlbound = computebalbounds(log, all_data, bus)

    # net sum of flows on branches = generation + stoc generation - correction  - demand
    # qtty on the right is what we represent using the variable InjPvar
    # Pubound and Plbound are the upper and lower bounds on generation - demand
    # we next adjust them to account for stoc generation

    stocvalue = 0
    if bus.nodeID in stocs.keys():
      stocvalue = stocs[bus.nodeID].average
    Plbound += stocvalue
    Pubound += stocvalue
    if stocvalue != 0:
      log.joint(" >> adjusted balance at " + str(bus.nodeID) + " by " + str(stocvalue) + " to account for stocs\n")

    if bus.nodetype == 3:
      Pubound = GRB.INFINITY
      Plbound = -GRB.INFINITY
    else:
      sumlower += Plbound

    
    InjPvar[bus] = themodel.addVar(obj = 0.0, lb = Plbound, ub = Pubound, name = "IP_"+str(bus.nodeID))
    
  log.joint("balance variables\n")
  
  for bus in buses.values():

    for genid in bus.genidsbycount:
      gen = gens[genid]
      lower = gen.Pmin*gen.status
      upper = gen.Pmax*gen.status
  
      if bus.nodetype == 3:
        upper = GRB.INFINITY
        lower = -GRB.INFINITY
      GenPvar[gen] = themodel.addVar(obj = 0.0, lb = lower, ub = upper, name = "GP_"+str(gen.count)+"_"+str(gen.nodeID))
      
  log.joint("generator variables\n")

  for branch in branches.values():
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]
    #Pvar
    ubound = branch.limit
    lbound = -branch.limit
    if all_data['flowsfilename'] != 'NONE' and branch.count in all_data['inputPf']:
      ubound = all_data['inputPf'][branch.count] + epsilon3
      lbound = all_data['inputPf'][branch.count] - epsilon3

    Pvar[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "P_"+str(branch.count)+","+str(f) + "," + str(t))
  log.joint("branch variables\n")

  lincostvar = themodel.addVar(obj = 1.0, lb = -GRB.INFINITY, ub = GRB.INFINITY, name = "lincost")
  quadcostvar = themodel.addVar(obj = 1.0, lb = 0, ub = GRB.INFINITY, name = "quadcost")

  constobjval = 0
  for gen in gens.values():  
    if gen.status > 0:
      constobjval += gen.costvector[2]
  
  # if the bus also participates in correction we adjust the constant objective alue
  corrvalue = 0
  S2 = all_data['sumvarstocs']
  log.joint(" sumvar of stocs: " + str(S2) + "\n")
  for bus in buses.values():
    for genid in bus.genidsbycount:
      if genid in alphas.keys():
        alphavalue = all_data['alphas'][genid].fixedvalue
        thisvalue = S2*alphavalue*alphavalue*gens[genid].costvector[0]
        corrvalue += thisvalue
        log.joint(">>>>>>> at buscount " + str(bus.count) + " generator id " + str(genid) + ": ")
        log.joint("  alphasq " + str(alphavalue*alphavalue) + ", ")
        log.joint("  costvector0 " + str(gens[genid].costvector[0]) + "\n")
        log.joint(">>>>>>> so correction term " + str(thisvalue) + "\n")

        #corrvalue -= alphas[genid].fixedvalue*all_data['sumavestocs']
        #when alpha is not fixed this expression ^^ has to be added to the constraint

  constobjval += corrvalue
  log.joint(" constant in objective: " + str(constobjval) + "\n")
    
  
  constvar = themodel.addVar(obj = constobjval, lb = 1.0, ub = 1.0, name = "constant")

  themodel.update()

  #constraints

  #cost definitions
  
  qcost = QuadExpr()
  
  for gen in gens.values():
    if gen.costvector[0] != 0:
      qcost += gen.costvector[0]*GenPvar[gen]*GenPvar[gen]
  
  themodel.addConstr(qcost <= quadcostvar, name = "qcostdef")

  coeff = [gen.costvector[1] for gen in gens.values()]
  variables = [GenPvar[gen] for gen in gens.values()]
  expr = LinExpr(coeff, variables)
  themodel.addConstr(expr == lincostvar, name="lincostdef")

  log.joint("cost def constraints\n")

  #define flow variables

  for branch in branches.values():
    f = branch.f
    t = branch.t
    coeff = 1/(branch.x*branch.ratio)
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]

    constrname = "Pdef_"+str(branch.count)+","+str(f)+","+str(t)
    expr = LinExpr()
    expr += coeff*thetavar[buses[count_of_f]]
    expr += -coeff*thetavar[buses[count_of_t]]
    themodel.addConstr(Pvar[branch] == expr - coeff*branch.angle_rad, name = constrname)

  log.joint("flow var def constraints\n")

  
  for bus in buses.values():
    nodeID = bus.nodeID
    constrname = "PBaldef"+str(nodeID)
    outflowexpr = LinExpr()
    for branchid in bus.frombranchids.values():
      outflowexpr += Pvar[ branches[branchid] ]
    for branchid in bus.tobranchids.values():
      outflowexpr -= Pvar[ branches[branchid] ]
      
    themodel.addConstr(outflowexpr == InjPvar[bus], name = constrname)

  log.joint("balance constraints\n")

  for bus in buses.values():

    stocvalue = 0
    if bus.nodeID in stocs.keys():
      stocvalue = stocs[bus.nodeID].average

    if len(bus.genidsbycount) > 0:
      constrname = "Bus_PInjdef_"+str(bus.nodeID)
      expr = LinExpr()
      for genid in bus.genidsbycount:
        gen = gens[genid]
        expr -= GenPvar[gen]
      expr += InjPvar[bus]
      themodel.addConstr(expr == stocvalue - bus.Pd, name = constrname)

  log.joint("generator def constraints\n")      

  themodel.update()

  log.joint("writing to lpfile " + all_data['lpfilename'] + "\n")
  themodel.write(all_data['lpfilename'])
  themodel.optimize()

  if themodel.status == GRB.status.INF_OR_UNBD:
    log.joint('->LP infeasible or unbounded\n')
    return themodel.status, 0

  if themodel.status == GRB.status.INFEASIBLE:
    log.joint('->LP infeasible\n')
    return themodel.status, 0

  if themodel.status == GRB.status.UNBOUNDED:
    log.joint('->LP unbounded\n')
    return themodel.status, 0

  log.joint(' objective = %g\n' % themodel.objVal)
  loud = 0
  totalqcost = totallcost = 0
  for bus in buses.values():
    if len(bus.genidsbycount) > 0:
      for genid in bus.genidsbycount:
        gen = gens[genid]
        genvalue = GenPvar[gen].x
        lincosterm = gen.costvector[1]*genvalue
        quadcosterm = gen.costvector[0]*genvalue*genvalue
        costerm = lincosterm + quadcosterm
        totalqcost += quadcosterm
        totallcost += lincosterm
        if loud:
          log.joint(" bus " + str(bus.nodeID) + " gen " + str(genid) + " -> " + str(genvalue) + " qcost " + str(quadcosterm) + " lcost " + str(lincosterm) + "\n")

  log.joint(" total qcost: " + str(totalqcost) + "\n")
  log.joint(" total lcost: " + str(totallcost) + "\n")
  totalvariablecost = totalqcost + totallcost
  log.joint(" total variable cost: " + str(totalvariablecost) + "\n")

  
  loud = 0

  psol = {}
  for geni in gens.values():
    pi = GenPvar[geni].x
    psol[geni] = pi

  costestimator = 0
  for geni in gens.values():
    i = geni.count
    if geni.count in all_data['alphas']:
      alphai = all_data['alphas'][geni.count].fixedvalue
      c0i = geni.costvector[0]
      c1i = geni.costvector[1]
      pi = psol[geni]

      costestimator += c0i*(pi*pi + alphai*alphai*S2) + c1i*pi

  log.joint(" costestimator = " + str(costestimator) + "\n")


  
  
    
def lpformulator_ac(log, all_data):
  log.joint("AC formulation\n")

  starttime = time.time()

  buses = all_data['buses']
  branches = all_data['branches']
  gens = all_data['gens']
  IDtoCountmap = all_data['IDtoCountmap']

  themodel = Model("csmodel")

  cvar = {}
  svar = {}
  Pvar_f = {}
  Qvar_f = {}
  Pvar_t = {}
  Qvar_t = {}
  BalPvar = {}
  BalQvar = {}
  GenPvar = {}
  GenQvar = {}
  
  epsilon_3 = all_data['epsilon_3']

  
  if all_data['voltsfilename'] != 'NONE':
    reader.generateinputcs(log,all_data)
    reader.generateinputeandf(log,all_data)
      
  #generate c variables
  for bus in buses.values():
    maxprod = bus.Vmax*bus.Vmax
    minprod = bus.Vmin*bus.Vmin
          
    ubound = maxprod
    lbound = minprod
          
    if all_data['voltsfilename'] != 'NONE' and bus.nodeID in all_data['inputcc']:
      ubound = all_data['inputcc'][bus.nodeID] + epsilon_3
      lbound = all_data['inputcc'][bus.nodeID] - epsilon_3

    cvar[bus] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "c_"+str(bus.nodeID)+","+str(bus.nodeID))
    Pubound, Plbound, Qubound, Qlbound = computebalbounds(log, all_data, bus)

    
    BalPvar[bus] = themodel.addVar(obj = 0.0, lb = Plbound, ub = Pubound, name = "IP_"+str(bus.nodeID))
    BalQvar[bus] = themodel.addVar(obj = 0.0, lb = Qlbound, ub = Qubound, name = "IQ_"+str(bus.nodeID))

    

    for genid in bus.genidsbycount:
      gen = gens[genid]
      lower = gen.Pmin*gen.status
      upper = gen.Pmax*gen.status
      if bus.nodetype == 3:
        upper = GRB.INFINITY
        lower = -GRB.INFINITY
      GenPvar[gen] = themodel.addVar(obj = 0.0, lb = lower, ub = upper, name = "GP_"+str(gen.count)+"_"+str(gen.nodeID))
      lower = gen.Qmin*gen.status
      upper = gen.Qmax*gen.status
      if bus.nodetype == 3:
        upper = GRB.INFINITY
        lower = -GRB.INFINITY
      GenQvar[gen] = themodel.addVar(obj = 0.0, lb = lower, ub = upper, name = "GQ_"+str(gen.count)+"_"+str(gen.nodeID))
    

  
  for branch in branches.values():
      f = branch.f
      t = branch.t
      count_of_f = IDtoCountmap[f]
      count_of_t = IDtoCountmap[t]
      maxprod = buses[count_of_f].Vmax*buses[count_of_t].Vmax
      minprod = buses[count_of_f].Vmin*buses[count_of_t].Vmin

      #first, c
      ubound = ubasic = maxprod
      lbound = lbasic = minprod*math.cos(branch.maxangle_rad)
      #stupid matpower

      if branch.upperanglenone == 1:
        ubound = maxprod
        lbound = 0

      if branch.loweranglenone == 1:
        ubound = maxprod
        lbound = 0

      
      if all_data['voltsfilename'] != 'NONE' and (f,t) in all_data['inputcs']:
          ubound = all_data['inputcs'][f,t][0] + epsilon_3
          lbound = all_data['inputcs'][f,t][0] - epsilon_3
          if lbound > ubasic or ubound < lbasic:
              log.joint("forced bounds for branch "+str(branch.count)+" "+str(f) + " " + str(t) + " are illegal\n")
              log.joint(" per forced_u " + str(ubound) + " forced_l " + str(lbound) + "\n")
              log.joint("  basic: " + str(ubasic) + " " + str(lbasic) + "\n")
              log.joint("  maxprod " + str(maxprod) + " minprod " + str(minprod) + " cosmaxangle "  + str(math.cos(branch.maxangle_rad)) + "\n")
              log.stateandquit(" quitting\n")
                                                                         
                                          
      cvar[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "c_"+str(branch.count)+","+str(f) + "," + str(t))

      #next, s

      if branch.maxangle_rad > 0:      
        ubound = maxprod*math.sin(branch.maxangle_rad)
      else:
        ubound = minprod*math.sin(branch.maxangle_rad)        
      if branch.minangle_rad <= 0:
        lbound = maxprod*math.sin(branch.minangle_rad)
      else:
        lbound = minprod*math.sin(branch.minangle_rad)


      if branch.upperanglenone == 1:
        ubound = maxprod
      if branch.loweranglenone == 1:
        lbound = -maxprod

      # print "f", branch.f, "t", branch.t, "ubound", ubound, "lbound", lbound, "none", branch.upperanglenone

      # breakexit("poo")


      if all_data['voltsfilename'] != 'NONE' and (f,t) in all_data['inputcs']:
          ubound = all_data['inputcs'][f,t][1]  + epsilon_3
          lbound = all_data['inputcs'][f,t][1] - epsilon_3
                                  
      svar[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "s_"+str(branch.count)+","+str(f) + "," + str(t))


      
  themodel.update()

  for branch in branches.values():
      f = branch.f
      t = branch.t
      count_of_f = IDtoCountmap[f]
      count_of_t = IDtoCountmap[t]

      #P and Q both ways  

      ubound = branch.limit
      lbound = -branch.limit

      if all_data['flowsfilename'] != 'NONE' and branch.count in all_data['inputPf']:
        ubound = all_data['inputPf'][branch.count] + epsilon_3
        lbound = all_data['inputPf'][branch.count] - epsilon_3

      Pvar_f[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "P_"+str(branch.count)+","+str(f) + "," + str(t))
      if all_data['flowsfilename'] != 'NONE' and branch.count in all_data['inputPt']:
        ubound = all_data['inputPt'][branch.count] + epsilon_3
        lbound = all_data['inputPt'][branch.count] - epsilon_3
      
      Pvar_t[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "P_"+str(branch.count)+","+str(t) + "," + str(f))

      ubound = branch.limit
      lbound = -branch.limit


      if all_data['flowsfilename'] != 'NONE' and branch.count in all_data['inputQf']:
        ubound = all_data['inputQf'][branch.count] + epsilon_3
        lbound = all_data['inputQf'][branch.count] - epsilon_3

      Qvar_f[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "Q_"+str(branch.count)+","+str(f) + "," + str(t))

      if all_data['flowsfilename'] != 'NONE' and branch.count in all_data['inputQt']:
        ubound = all_data['inputQt'][branch.count] + epsilon_3
        lbound = all_data['inputQt'][branch.count] - epsilon_3
      
      Qvar_t[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "Q_"+str(branch.count)+","+str(t) + "," + str(f))

  lincostvar = themodel.addVar(obj = 1.0, lb = -GRB.INFINITY, ub = GRB.INFINITY, name = "lincost")
  quadcostvar = themodel.addVar(obj = 1.0, lb = 0, ub = GRB.INFINITY, name = "quadcost")

  constobjval = 0  
  for gen in gens.values():  
    if gen.status > 0:
      constobjval += gen.costvector[ gen.costdegree ]
  
  constvar = themodel.addVar(obj = constobjval, lb = 1.0, ub = 1.0, name = "constant")

  themodel.update()

  #constraints

  #lincost definition
  
  coeff = [gen.costvector[gen.costdegree-1] for gen in gens.values()]
  variables = [GenPvar[gen] for gen in gens.values()]
  expr = LinExpr(coeff, variables)
  themodel.addConstr(expr == lincostvar, name="lincostdef")

  qcost = QuadExpr()
  for gen in gens.values():
    if gen.costdegree == 2 and gen.costvector[0] != 0:
      qcost += gen.costvector[0]*GenPvar[gen]*GenPvar[gen]
  
  themodel.addConstr(qcost <= quadcostvar, name = "qcostdef")

  #define flow variables

  for branch in branches.values():
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]

    #  Gff cff + Gft cft + Bft sft
    constrname = "Pdef_"+str(branch.count)+","+str(f)+","+str(t)
    expr = LinExpr()
    expr += branch.Gff*cvar[buses[count_of_f]]
    expr += branch.Gft*cvar[branch]
    expr += branch.Bft*svar[branch]
    themodel.addConstr(expr == Pvar_f[branch], name = constrname)

    #  Gtt ctt + Gtf cft + Btf stf = Gtt ctt + Gtf cft - Btf sft
    constrname = "Pdef_"+str(branch.count)+","+str(t)+","+str(f)
    expr = LinExpr()
    expr += branch.Gtt*cvar[buses[count_of_t]]
    expr += branch.Gtf*cvar[branch]
    expr += -branch.Btf*svar[branch] #minus because svarft = -svartf
    themodel.addConstr(expr == Pvar_t[branch], name = constrname)

  for branch in branches.values():
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]
    constrname = "Qdef_"+str(branch.count)+","+str(f)+","+str(t)

    # -Bff cff - Bft cft + Gft sft
    expr = LinExpr()
    expr += -branch.Bff*cvar[buses[count_of_f]]
    expr += -branch.Bft*cvar[branch]
    expr += +branch.Gft*svar[branch]
    themodel.addConstr(expr == Qvar_f[branch], name = constrname)

    #  -Btt ctt - Btf cft + Gtf stf = -Btt ctt - Btf cft - Gtf sft 
    constrname = "Qdef_"+str(branch.count)+","+str(t)+","+str(f)
    expr = LinExpr()
    expr += -branch.Btt*cvar[buses[count_of_t]]
    expr += -branch.Btf*cvar[branch]
    expr += -branch.Gtf*svar[branch] #again, same minus
    themodel.addConstr(expr == Qvar_t[branch], name = constrname)

 
  for bus in buses.values():
    constrname = "PBaldef"+str(bus.nodeID)
    expr = LinExpr()
    for branchid in bus.frombranchids.values():
      expr += Pvar_f[ branches[branchid] ]
    for branchid in bus.tobranchids.values():
      expr += Pvar_t[ branches[branchid] ]


    if bus.Gs != 0:
      expr += bus.Gs*cvar[bus]
      
      
    themodel.addConstr(expr == BalPvar[bus], name = constrname)

  for bus in buses.values():
    constrname = "QBaldef"+str(bus.nodeID)
    expr = LinExpr()
    for branchid in bus.frombranchids.values():
      expr += Qvar_f[ branches[branchid] ]
    for branchid in bus.tobranchids.values():
      expr += Qvar_t[ branches[branchid] ]
 
    if bus.Bs != 0:
      expr += (-bus.Bs)*cvar[bus]

    themodel.addConstr(expr == BalQvar[bus], name = constrname)

  for branch in branches.values():
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]
    trigexp = QuadExpr()
    constrname = "quad_"+str(branch.count)+"_"+str(f)+","+str(t)

    trigexp += cvar[branch]*cvar[branch] + svar[branch]*svar[branch] - cvar[buses[count_of_f]]*cvar[buses[count_of_t]]

    themodel.addConstr(trigexp <= 0, name = constrname)

  for branch in branches.values():
    if branch.constrainedflow:
      f = branch.f
      t = branch.t
      constrname = "limit_f_"+str(branch.count)+","+str(f)+","+str(t)
      limexp = QuadExpr()
      limexp += Pvar_f[branch]*Pvar_f[branch] + Qvar_f[branch]*Qvar_f[branch]
      themodel.addConstr(limexp <= branch.limit**2, name = constrname)

      
      constrname = "limit_t_"+str(branch.count)+","+str(t)+","+str(f)
      limexp = QuadExpr()
      limexp += Pvar_t[branch]*Pvar_t[branch] + Qvar_t[branch]*Qvar_t[branch]
      #themodel.cbLazy(limexp <= branch.limit**2)
      themodel.addConstr(limexp <= branch.limit**2, name = constrname)

  for bus in buses.values():

    if len(bus.genidsbycount) > 0:
      constrname = "Bus_PInj_"+str(bus.nodeID)
      expr = LinExpr()
      for genid in bus.genidsbycount:
        gen = gens[genid]
        expr += GenPvar[gen]
      expr += -BalPvar[bus]
      themodel.addConstr(expr == bus.Pd, name = constrname)

      constrname = "Bus_QInj_"+str(bus.nodeID)
      expr = LinExpr()
      for genid in bus.genidsbycount:
        gen = gens[genid]
        expr += GenQvar[gen]
      expr += -BalQvar[bus]
      themodel.addConstr(expr == bus.Qd, name = constrname)

  pi = math.pi
  for branch in branches.values():
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]

    if branch.maxangle_rad < pi/2 and branch.upperanglenone == 0:
      expr = LinExpr()
      constrname = "maxangle_"+str(branch.count)+"_"+str(f)+","+str(t)
      expr += svar[branch] - math.tan(branch.maxangle_rad)*cvar[branch]
      themodel.addConstr(expr <= 0, name = constrname)

      #print "count", branch.count, "f", branch.f, "t", branch.t, "maxangle", branch.maxangle_rad


    if branch.minangle_rad > -pi/2  and branch.loweranglenone == 0:
      expr = LinExpr()
      constrname = "minangle_"+str(branch.count)+"_"+str(f)+","+str(t)
      expr += svar[branch] - math.tan(branch.minangle_rad)*cvar[branch]
      themodel.addConstr(expr >= 0, name = constrname)
      #print f, t, branch.minangle, branch.minangle_rad, math.tan(branch.minangle_rad)

  if all_data['budgetrow'] == 'YES':
    #add budget row
    expr = LinExpr()
    constrname = "budgetrow"
    expr += lincostvar + quadcostvar + constobjval*constvar
    objexpr = expr
    themodel.addConstr(expr <= all_data['budget'], name = constrname)      
      
  themodel.update()

  endtime = time.time()

  log.joint(' formulation time: %g\n' %(endtime - starttime))

  log.joint("writing to lpfile " + all_data['lpfilename'] + "\n")  
  themodel.write(all_data['lpfilename'])
  log.joint("writing to mpsfile\n")
  themodel.write('foo.mps')


  all_data['themodel'] = themodel

  #get vars in there
  all_data['cvar'] = cvar
  all_data['svar'] = svar

  if all_data['digitize']:
    digitizer(log, all_data)

  themodel.Params.BarHomogeneous = 1

  if all_data['obbt'] == 'YES':

    code, oldvalue = gur_optimize(log, themodel, 0)
    log.joint("opt code " + str(code) + " value " + str(oldvalue) + "\n")

    #budgconst = themodel.getConstrByName("budgetrow")
    #budgconst.rhs = oldvalue
    breakexit("preopt")


    itcount = 0
    maxitcount = 10
    gap = 1e20

    themodel.write('obbt_'+str(itcount) + '.mps')

    itcount += 1
    

    #while gap > epsilon4 and itcount < maxitcount:
    while itcount < maxitcount:      
    
      obbt(log, all_data, themodel, itcount)

      #now resolve
      themodel.setObjective(objexpr, GRB.MINIMIZE)
      themodel.write('obbt_'+str(itcount) + '.mps')
      code, value = gur_optimize(log, themodel, 0)
      gap = value - oldvalue
      log.joint("it " + str(itcount) + " reopt code " + str(code) + " value " + str(value) + " gap " + str(gap) +"\n")

      oldvalue = value

      itcount +=1

    
    breakexit("done with obbt")
  themodel.optimize()

  if themodel.status == GRB.status.INF_OR_UNBD:
    log.joint('->LP infeasible or unbounded\n')
    return themodel.status, 0

  if themodel.status == GRB.status.INFEASIBLE:
    log.joint('->LP infeasible\n')
    return themodel.status, 0

  if themodel.status == GRB.status.UNBOUNDED:
    log.joint('->LP unbounded\n')
    return themodel.status, 0

  log.joint(' objective = %g\n' % themodel.objVal)


  if all_data['doanalysis'] == 0:
    return themodel.status, themodel.objVal

  cvalues = {}
  svalues = {}

  for bus in buses.values():
    log.joint( str(cvar[bus].varname) + " "  +str(cvar[bus].x) + "\n")
    cvalues[bus] = cvar[bus].x

  for branch in branches.values():
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]
    log.joint( str(cvar[branch].varname) + " " + str(cvar[branch].x) + "\n")
    cvalues[branch] = cvar[branch].x
    log.joint( str(svar[branch].varname) + " " + str(svar[branch].x) + "\n")
    svalues[branch] = svar[branch].x
    sin_estimate = svalues[branch]/(cvalues[buses[count_of_f]]*cvalues[buses[count_of_t]])**.5
    rangle_estimate = math.asin(sin_estimate)
    angle_estimate = rangle_estimate*180/math.pi
    log.joint(" from " + str(f) + " to " + str(t) +", sin_estimate: " + str(sin_estimate) + " angle_estimate: " + str(angle_estimate) + "\n")

  count = 0

  hit = {}
  for branch in branches.values():
    lhs = cvalues[branch]**2 + svalues[branch]**2
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]
    #print lhs, count_of_f, count_of_t

    rhs = cvalues[buses[count_of_f]]*cvalues[buses[count_of_t]]

    log.joint("branch " + str(branch.count) + ", " + str(f) + " " + str(t) + ", lhs " +str(lhs) + " rhs "  +str(rhs) + " gap  " + str((rhs - lhs)/rhs) + "\n")
    if (rhs - lhs)/rhs > 1.0e-04:
      count += 1
      hit[f] = hit[t] = 1
      log.joint(" hit " + str(f) + " and " + str(t) + "\n")


  hitcount = 0
  for bus in buses.values():
    if bus.nodeID in hit:
      log.joint(" bus " + str(bus.nodeID) + " hit\n")
      hitcount += 1
  log.joint(" inactive buses " +str(count) + " branches "  + str(len(branches)) + " buses "  + str(len(buses)) + " buses hit by slack branches "  + str(hitcount) + "\n")

  return themodel.status, themodel.objVal

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


def obbt(log, all_data, themodel, itcount):
  #breakexit("in obbt")

  # change objective coefficient of lincost, quadcost and constant to zero
  themodel.setParam("OutputFlag", False)


  numhi = numlo = 0

  minwidth = 1e20
  badvar = 'none'
  
  for var in themodel.getVars():

    if var.varname == 'lincost' or var.varname == 'quadcost' or var.varname == 'constant':
      continue
    
    themodel.setObjective(var, GRB.MAXIMIZE)
    #themodel.write('foo.lp')
    #breakexit("huh")
    code, value = gur_optimize(log, themodel, 1)
    themodel.write('obbt'+'_'+var.varname+'_'+'up'+str(itcount)+'.mps')
    if code == 0:
      if value < var.ub - epsilon4:
        log.joint("variable " + var.varname)
        log.joint(" value " + str(value))
        log.joint(" old ub was " + str(var.ub) + "\n")
        var.ub = value
        width = var.ub - var.lb
        if width < minwidth:
          minwidth = width
          badvar = var.varname
        numhi += 1
        themodel.update()
    else:
      log.joint(' bad up opt code for ' + var.varname + ' lb ' + str(var.lb) + ' ub ' + str(var.ub) + '\n')
      log.joint(' minwidth ' + str(minwidth) + ' badvar ' + badvar + '\n')

    # next do the min
    # change objective sense
    themodel.setObjective(var, GRB.MINIMIZE)
    #themodel.write('foo.lp')
    #breakexit("huh")
    code, value = gur_optimize(log, themodel, 1)
    themodel.write('obbt'+'_'+var.varname+'_'+'dn'+str(itcount)+'.mps')    
    if code == 0:

      if value > var.lb + epsilon4:
        log.joint("variable " + var.varname)
        log.joint(" value " + str(value))
        log.joint(" old lb was " + str(var.lb) + "\n")
        var.lb = value
        width = var.ub - var.lb
        if width < minwidth:
          minwidth = width
          badvar = var.varname
        numlo += 1
        themodel.update()
    else:
      log.joint(' bad dn opt code for ' + var.varname + ' lb ' + str(var.lb) + ' ub ' + str(var.ub) + '\n')
      log.joint(' minwidth ' + str(minwidth) + ' badvar ' + badvar + '\n')

  log.joint("numlo " + str(numlo) + " numhi " + str(numhi) + "\n")
  #breakexit("huh3")

  # restore objective coefficient of lincost, quadcost and constant 
  #breakexit("end of obbt")


def gur_optimize(log, themodel, nonleeway):

  
  themodel.params.method = 2

  themodel.optimize()


  if themodel.status == GRB.status.INF_OR_UNBD:
    log.joint('->LP infeasible or unbounded\n')
    return 1, 0
  elif themodel.status == GRB.status.INFEASIBLE:
    log.joint('->LP infeasible\n')
    return 1, 0
  elif themodel.status == GRB.status.UNBOUNDED:
    log.joint('->LP unbounded\n')
    return 1, 0
  elif themodel.status == GRB.status.OPTIMAL or nonleeway == 0:
    return 0, themodel.objVal
  else:
    log.joint(' --> other status ' + str(themodel.status) + '\n')
    return 1, 0

  return 0, themodel.objVal


  
