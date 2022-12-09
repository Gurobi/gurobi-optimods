import math
import cmath
import time

from myutils import *
from log import *

class Bus:
    def __init__(self, count, nodeID, nodetype, Pd, Qd, Gs, Bs, Vbase, Vmax,
                 Vmin, busline0):
        self.genidsbycount = []
        self.frombranchids = {}
        self.tobranchids   = {}
        self.count         = count
        self.nodeID        = nodeID
        self.nodetype      = nodetype
        self.Pd            = Pd
        self.Qd            = Qd
        self.Gs            = Gs
        self.Bs            = Bs
        self.Vbase         = Vbase
        self.Vmax          = Vmax
        self.Vmin          = Vmin
        self.busline0      = busline0
        self.inputvoltage  = False
        self.cffvarind     = -1
        self.Pinjvarind    = -1
        self.Qinjvarind    = -1
        self.Pbalance      = 0
        self.inputV        = 0
        self.outdegree = self.indegree = self.degree = 0

    def getbusline0(self):
        return self.busline0

    def addgenerator(self, log, generatorcount, generator):
        self.genidsbycount.append(generatorcount)
        #loud = 0
        #if loud:
        #  log.joint(" added generator # " + str(generatorcount) + " to bus ID " + str(self.nodeID))
        #  log.joint(" Pmax " + str(generator.Pmax) + " Pmin " + str(generator.Pmin) + "\n")

    def addfrombranch(self, log, id):
        quant = len(self.frombranchids)
        self.frombranchids[quant] = id
        self.outdegree += 1
        self.degree += 1

    def addtobranch(self, log, id):
        quant = len(self.tobranchids)
        self.tobranchids[quant] = id
        self.indegree += 1
        self.degree += 1

class Branch:
    def __init__(self, log, count, f, id_f, t, id_t, r, x, bc, rateAmva,
                 rateBmva, rateCmva, ratio, angle, maxangle, minangle,
                 status, defaultlimit, branchline0):
        self.count           = count
        self.f               = f
        self.t               = t
        self.id_f            = id_f
        self.id_t            = id_t
        self.r               = r
        self.x               = x
        self.bc              = bc
        self.count           = count
        self.branchline0     = branchline0
        self.rateAmva        = rateAmva
        self.rateBmva        = rateBmva
        self.limit           = rateAmva
        self.unboundedlimit  = False
        self.constrainedflow = 1

        if self.limit == 0:
            self.limit = defaultlimit
            self.constrainedflow = 0
            self.unboundedlimit = True

        if ratio == 0:
            ratio = 1

        self.ratio          = ratio
        self.angle          = angle
        self.angle_rad      = math.pi*angle/180.0
        self.maxangle       = maxangle
        self.maxangle_rad   = math.pi*maxangle/180.0
        self.minangle       = minangle
        self.minangle_rad   = math.pi*minangle/180.0
        self.upperanglenone = 0

        if maxangle == 360 or maxangle == 0:
            self.maxangle_rad   = 2*math.pi
            self.upperanglenone = 1

        self.loweranglenone = 0

        if minangle == -360 or minangle == 0:
            self.minangle_rad   = -2*math.pi
            self.loweranglenone = 1

        self.invratio2 = invratio2 = 1/ratio**2
        self.multtf    = multtf = 1/(ratio*cmath.exp(1j*self.angle_rad))
        self.multft    = multft = 1/(ratio*cmath.exp(-1j*self.angle_rad))
        #print 'multtf', multtf
        self.status = status
        self.z = z = r + x*1j
        self.y = y = 1/z
        self.Yff = Yff = (y + bc/2*1j)*invratio2
        self.Yft = Yft = -y*multft
        self.Ytf = Ytf = -y*multtf
        self.Ytt = Ytt = y + bc/2*1j
        self.Gff = Gff = (self.Yff).real
        self.Bff = Bff = (self.Yff).imag
        self.Gft = Gft = (self.Yft).real
        self.Bft = Bft = (self.Yft).imag
        self.Gtf = Gtf = (self.Ytf).real
        self.Btf = Btf = (self.Ytf).imag
        self.Gtt = Gtt = (self.Ytt).real
        self.Btt = Btt = (self.Ytt).imag

        self.inputcs   = False
        self.inputc    = 2
        self.inputs    = 2  #bogus #FIXME
        self.cftvarind = -1
        self.sftvarind = -1
        self.Pftvarind = -1
        self.Qftvarind = -1

#       loud = False
#       if self.angle_rad != 0:
#         loud = 1
#       if loud:
#        log.joint("\nbr " + str(count) + " f " + str(f) + " t " + str(t) +"\n")
#        log.joint("   idf " +  str(id_f) + " idt " + str(id_t) +"\n")
#        log.joint("   r " + str(r) + " x " + str(x) + " bb " + str(bc) +"\n")
#        log.joint("   ratio " + str(self.ratio) + " angle " + str(angle) + " angle_rad: " + str(self.angle_rad) + "\n")
#        log.joint("   y " + str(y) + "\n")
#        log.joint("       Yff " + str(Yff) + " , Yft " + str(Yft) + " , Ytf " + str(Ytf) +  " , Ytt " + str(Ytt) + "\n")

    def getbranchline0(self):
        return self.branchline0

    def show(self, log):
        log.joint(" < " + str(self.f) + " , " + str(self.t) + " > ")
        log.joint(" r " + str(self.r) + " x " + str(self.x) + " bc " + str(self.bc))
        log.joint("\n")
        log.joint(" ra " + str(self.ratio) + " ang " + str(self.angle) )
        log.joint("\n")

class Gen:
    def __init__(self, count, nodeID, Pg, Qg, status, Pmax, Pmin, Qmax,
                 Qmin, line0):
        self.count       = count
        self.nodeID      = nodeID
        self.Pg          = Pg
        self.Qg          = Qg
        self.status      = status
        self.Pmax        = Pmax
        self.Pmin        = Pmin
        self.Qmax        = Qmax
        self.Qmin        = Qmin
        self.line0       = line0
        self.costlinenum = -1
        self.Pvarind     = -1
        self.Qvarind     = -1

    def addcost(self, log, costvector, linenum):
        self.costvector  = costvector
        self.costdegree  = len(costvector) - 1
        self.costlinenum = linenum

    def getline0(self):
        return self.line0

def readcase(alldata):
    log          = alldata['log']
    casefilename = alldata['casefilename']
    log.joint("Reading case file %s\n"%casefilename)

    starttime = time.time()

    try:
        f     = open(casefilename, "r")
        lines = f.readlines()
        f.close()
    except:
        log.stateandquit("Error: Cannot open file %s\n"%casefilename)

    readcase_thrulines(log, alldata, lines)
    alldata['casefilelines'] = lines

    endtime = time.time()

    log.joint("Reading time: %f s\n"%(endtime - starttime))

def readcase_thrulines(log, alldata, lines):
    log.joint("Reading case file thru lines\n")

    lookingforbus     = 1
    linenum           = 2
    numlines          = len(lines)
    numisolated       = 0
    alldata['refbus'] = -1

    while linenum <= numlines:
        line     = lines[linenum-1]
        thisline = line.split()

        if len(thisline) > 0:
          theword = thisline[0]
          if theword[0:4] == 'mpc.':
            log.joint("Found %s on line %d\n"%(theword, linenum))
            if theword == "mpc.baseMVA":
              boo = thisline[2]
              if boo[len(boo)-1] == ";":
                boo = boo[:len(boo)-1]
              alldata['baseMVA'] = float(boo)
              baseMVA = alldata['baseMVA']
              log.joint(" baseMVA: %s\n"%str(baseMVA))
            elif theword == 'mpc.bus':
              buses              = {}
              IDtoCountmap       = {}
              lookingforendofbus = 1
              numbuses           = 0
              slackbus           = -1
              linenum += 1

              numPload = 0
              sumload  = 0
              sumPd = sumQd = 0
              while lookingforendofbus and linenum <= numlines:
                line     = lines[linenum-1]
                thisline = line.split()
                length   = len(thisline)

                if thisline[0] == "];":
                  log.joint(" Found end of bus section on line %d\n"%linenum)
                  lookingforendofbus = 0
                  break

                numbuses += 1
                if thisline[1] == "3":
                  slackbus = int(thisline[0])
                  log.joint(" Slack bus: %d"%slackbus)

                if thisline[0] != '%':
                  nodeID, nodetype = int(thisline[0]), int(thisline[1])

                  if nodetype != 1 and nodetype != 2 and nodetype != 3 and nodetype != 4:
                    log.stateandquit("Error: Bad bus %s has type %s\n"%(thisline[0], thisline[1]))

                  if nodetype == 4:
                    #log.joint("bus " + thisline[0] + " is isolated\n")
                    numisolated += 1

                  # trim ; if present
                  Vmin = thisline[12]
                  if Vmin[len(Vmin)-1] == ';':
                    Vmin = Vmin[:len(Vmin)-1]

                  Vmin  = float(Vmin)
                  Pd    = float(thisline[2])
                  Qd    = float(thisline[3])
                  Gs    = float(thisline[4])
                  Bs    = float(thisline[5])
                  Vbase = float(thisline[9])
                  Vmax  = float(thisline[11])

                  buses[numbuses] = Bus(numbuses, nodeID, nodetype, Pd/baseMVA,
                                        Qd/baseMVA, Gs/baseMVA, Bs/baseMVA,
                                        Vbase, Vmax, Vmin, linenum-1)

                  if nodetype == 3:
                    log.joint('Bus %d ID %d is the reference bus\n'%(numbuses, nodeID))
                    alldata['refbus'] = numbuses
                    #breakexit('ref')

                  if nodetype == 1 or nodetype == 2 or nodetype == 3:
                    sumPd += Pd
                    sumQd += Qd

                  IDtoCountmap[nodeID] = numbuses
                  numPload += (Pd > 0)

                linenum += 1

              alldata['buses']        = buses
              alldata['numbuses']     = numbuses
              alldata['sumPd']        = sumPd
              alldata['sumQd']        = sumQd
              alldata['IDtoCountmap'] = IDtoCountmap
              alldata['slackbus']     = slackbus

              log.joint(" sumloadPd %f numPload %f\n"%(sumPd, numPload))
              log.joint(" sumloadQd %f\n"%sumQd)

              if lookingforendofbus:
                log.stateandquit("Error: Could not find bus data section\n")

              if slackbus < 0:
                log.joint(" Could not find slack bus\n")

              log.joint(" %d buses\n"%numbuses)
              if numisolated > 0:
                log.joint(" isolated buses: %d\n"%numisolated)

            elif theword == 'mpc.gen':
              gens               = {}
              lookingforendofgen = 1
              gencount           = 0
              linenum += 1

              summaxgenP = summaxgenQ = 0
              while lookingforendofgen and linenum <= numlines:
                line     = lines[linenum-1]
                thisline = line.split()

                if thisline[0] == "];":
                  alldata['endofgen'] = linenum
                  log.joint(" Found end of gen section on line %d\n"%linenum)
                  lookingforendofgen = 0
                  break

                gencount += 1

                nodeID = int(thisline[0])
                Pg     = float(thisline[1])
                Qg     = float(thisline[2])
                status = int(thisline[7])
                Pmax   = float(thisline[8])
                Pmin   = float(thisline[9])
                Qmax   = float(thisline[3])
                Qmin   = float(thisline[4])

                if status <= 0:
                  status = 0
                else:
                  status = 1

                #log.joint("generator in bus ID " + str(nodeID) )
                if nodeID in IDtoCountmap.keys():
                  idgencount = IDtoCountmap[nodeID]
                  gens[gencount] = Gen(gencount, nodeID, Pg, Qg, status,
                                       Pmax/baseMVA, Pmin/baseMVA, Qmax/baseMVA,
                                       Qmin/baseMVA, linenum-1)
                  buses[idgencount].addgenerator(log, gencount, gens[gencount])

                  if buses[idgencount].nodetype == 2 or buses[idgencount].nodetype == 3:  #but not 4
                    summaxgenP += Pmax
                    summaxgenQ += Qmax

                else:
                  log.stateandquit("Error: Generator # %d in noncexistent bus ID %d\n"%(gencount, nodeID))

                linenum += 1

              if lookingforendofgen:
                log.stateandquit("Error: Could not find end of generator section\n")

              alldata['gens']    = gens
              alldata['numgens'] = len(gens)
              busgencount        = 0
              for bus in buses.values():
                busgencount += len(bus.genidsbycount) > 0

              alldata['busgencount'] = busgencount
              log.joint(" number of generators: %d\n"%alldata['numgens'])
              log.joint(" number of buses with gens: %d\n"%alldata['busgencount'])
              alldata['summaxgenP'] = summaxgenP
              alldata['summaxgenQ'] = summaxgenQ              
              log.joint(" summaxPg: %f summaxQg: %f\n"%(summaxgenP, summaxgenQ))

            elif theword == 'mpc.branch':
              branches              = {}
              defaultlimit          = 1e20
              lookingforendofbranch = 1
              numbranches           = 0
              activebranches        = 0
              zerolimit             = 0
              linenum += 1

              while lookingforendofbranch and linenum <= numlines:
                line     = lines[linenum-1]
                thisline = line.split()

                if thisline[0] == "];":
                  log.joint(" Found end of branch section on line %d\n"%linenum)
                  lookingforendofbranch = 0
                  break

                numbranches += 1
                f        = int(thisline[0])
                t        = int(thisline[1])
                r        = float(thisline[2])
                x        = float(thisline[3])
                bc       = float(thisline[4])
                rateA    = float(thisline[5])
                rateB    = float(thisline[6])
                rateC    = float(thisline[7])
                ratio    = float(thisline[8])
                angle    = float(thisline[9])
                status   = int(thisline[10])
                minangle = float(thisline[11])
                maxangle = thisline[12]
                # trim ; at the end of line
                if maxangle[len(maxangle)-1] == ';':
                  maxangle = maxangle[:len(maxangle)-1]

                maxangle = float(maxangle)

                if maxangle < minangle:
                  log.stateandquit("Error: Branch # %d has illegal angle constraints\n"%numbranches)

                id_f = IDtoCountmap[f]
                id_t = IDtoCountmap[t]

                if status:
                  branches[numbranches] = Branch(log,numbranches, f, id_f, t,
                                                 id_t, r, x, bc, rateA/baseMVA,
                                                 rateB/baseMVA, rateC/baseMVA,
                                                 ratio, angle, maxangle, minangle,
                                                 status, defaultlimit, linenum-1)
                  zerolimit += (branches[numbranches].constrainedflow == 0)
                  activebranches += 1
                  buses[id_f].addfrombranch(log, numbranches)
                  buses[id_t].addtobranch(log, numbranches)

                linenum += 1

              if lookingforendofbranch:
                log.stateandquit("Error: Could not find end of branch section\n")

              alldata['branches']    = branches
              alldata['numbranches'] = numbranches
              log.joint(" numbranches: %d active %d\n"%(numbranches, activebranches))
              log.joint(" %d unconstrained\n"%zerolimit)

            elif theword == 'mpc.gencost':
              lookingforendofgencost = 1
              gencostcount           = 1
              linenum += 1

              while lookingforendofgencost and linenum <= numlines:
                line     = lines[linenum-1]
                thisline = line.split()

                if thisline[0] == "];":
                  log.joint(" Found end of gencost section on line %d\n"%linenum)
                  alldata['endofgencost'] = linenum
                  lookingforendofgencost  = 0
                  break

                if gencostcount <= gencount:
                  costtype = int(thisline[0])
                  if costtype != 2:
                    log.stateandquit("Error: Cost of generator %d is not polynomial\n"%gencostcount)

                  degree = int(thisline[3]) - 1

                  if degree > 2 or degree < 0:
                    log.stateandquit("Error: Degree of cost function for generator %d is illegal\n"%gencostcount)

                  costvector = [0 for j in range(degree+1)]

                  for j in range(degree+1):
                    tmp = thisline[4+j]
                    if tmp[len(tmp)-1] == ";":
                      tmp = tmp[:len(tmp)-1]

                    costvector[j] = float(tmp)
                    costvector[j] *= (baseMVA)**(degree - j)
                    # print "gen", gencostcount, j, costvector[j]
                  #print costvector

                  gens[gencostcount].addcost(log, costvector, linenum)

                else:
                  log.stateandquit("Error: Read %d gen costs but only %d generators\n"%(gencostcount, gencount))

                gencostcount += 1
                linenum += 1

              if lookingforendofgencost:
                log.stateandquit("Error: Could not find end of gencost section\n")

              linenum += 1

        linenum += 1

def readvoltsfile(log, alldata):
    voltsfilename = alldata['voltsfilename']
    IDtoCountmap  = alldata['IDtoCountmap']
    buses         = alldata['buses']

    try:
        f     = open(voltsfilename, "r")
        lines = f.readlines()
        log.joint("reading volts file %s\n"%voltsfilename)
        f.close()
    except:
        log.stateandquit("Error: Cannot open file %s\n"%voltsfilename)

    inputvolts = {}
    numread = 0

    for linenum in range(len(lines)):
        thisline = lines[linenum].split()

        if len(thisline) > 0:
          if thisline[0] == 'bus':
            angle_rad         = float(thisline[5])*math.pi/180
            busid             = int(thisline[1])
            inputvolts[busid] = (float(thisline[3]), angle_rad)
            numread += 1

          elif thisline[0] == 'END':
            break
          else:
            log.stateandquit("Error: Illegal input %s\n"%thisline[0])

    log.joint("Read %d input voltages\n"%numread)
    alldata['inputvolts'] = inputvolts

def readflowsfile(log, alldata):
    flowsfilename = alldata['flowsfilename']
    IDtoCountmap  = alldata['IDtoCountmap']
    buses         = alldata['buses']
    try:
        f     = open(flowsfilename, "r")
        lines = f.readlines()
        log.joint(" Reading flows file %s\n"%flowsfilename)
        f.close()
    except:
        log.stateandquit("Cannot open flows file %s"%flowsfilename)

    inputPf = {}
    inputQf = {}
    inputPt = {}
    inputQt = {}
    numread = 0

    baseMVA = alldata['baseMVA']

    for linenum in range(len(lines)):
        thisline = lines[linenum].split()
        if len(thisline) > 0:

          if thisline[0] == 'branch':
            branchid          = int(thisline[1])
            inputPf[branchid] = float(thisline[7])/baseMVA
            inputPt[branchid] = float(thisline[9])/baseMVA
            inputQf[branchid] = float(thisline[11])/baseMVA
            inputQt[branchid] = float(thisline[13])/baseMVA

            numread += 1
          elif thisline[0] == 'END':
            break
          else:
            log.stateandquit("Error: Illegal input %s on line %s\n"%(thisline[0], thisline))

    log.joint("Read %d input flows\n"%numread)

    alldata['inputPf'] = inputPf
    alldata['inputPt'] = inputPt
    alldata['inputQf'] = inputQf
    alldata['inputQt'] = inputQt

def writegv(log, alldata, gvfilename):
  try:
    f = open(gvfilename, "w")
    log.joint("Writing to gv file %s\n"%gvfilename)
  except:
    log.stateandquit("Cannot open file %s"%gvfilename)

  f.write("graph {\n")

  for bus in alldata['buses'].values():
    f.write("     " + str(bus.nodeID)+";\n")

  for branch in alldata['branches'].values():
    f.write("     " + str(branch.f)+" -- " + str(branch.t)+";\n")

  f.write("}\n")
  f.close()

def generateinputcs(log, alldata):
  log.joint("  generating input c,s values\n")

  inputcc      = {}
  inputcs      = {}
  inputvolts   = alldata['inputvolts']
  buses        = alldata['buses']
  branches     = alldata['branches']
  IDtoCountmap = alldata['IDtoCountmap']

  for busid in inputvolts:
    M              = inputvolts[busid][0]
    inputcc[busid] = M*M
    #log.joint(str(busid) +" is in input volts with M " +  str(M) + "\n")

  for branch in branches.values():
    f          = branch.f
    t          = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]

    if f in inputvolts and t in inputvolts:
      Mf           = inputvolts[f][0]
      Mt           = inputvolts[t][0]
      af           = inputvolts[f][1]
      at           = inputvolts[t][1]
      angle        = af - at
      inputcs[f,t] = (Mf*Mt*math.cos(angle), Mf*Mt*math.sin(angle))
      inputcs[t,f] = (Mf*Mt*math.cos(angle), -Mf*Mt*math.sin(angle))

  alldata['inputcc'] = inputcc 
  alldata['inputcs'] = inputcs 

def generateinputeandf(log, alldata):
  log.joint("  generating input e,f values\n")

  inputve    = {}
  inputvf    = {}
  inputvolts = alldata['inputvolts']
  buses        = alldata['buses']
  IDtoCountmap = alldata['IDtoCountmap']

  for busid in inputvolts:
    M              = inputvolts[busid][0]
    A              = inputvolts[busid][1]
    inputve[busid] = M*math.cos(A)
    inputvf[busid] = M*math.sin(A)
    #log.joint(str(busid) +" is in input volts with M " +  str(M) + "\n")

  alldata['inputve'] = inputve
  alldata['inputvf'] = inputvf

def readdigits(log, alldata):
  Lfilename = alldata['Lfilename']
  log.joint("reading L file %s\n"%Lfilename)

  buses = alldata['buses']

  try:
    f     = open(Lfilename, "r")
    lines = f.readlines()
    f.close()
  except:
    log.stateandquit("Error: Cannot open file %s\n")

  L            = {}
  IDtoCountmap = alldata['IDtoCountmap'] 

  for bus in buses.values():
    L[bus] = 0

  for linenum in range(len(lines)):
    line = lines[linenum].split()

    if line[0] == 'default':
      lvalue = int(line[1])
      log.joint(" default L: %d\n"%lvalue)

      for bus in buses.values():
        L[bus] = lvalue

    elif line[0] == 'END':
      break
    else:
      ind = int(line[0])
      L[buses[IDtoCountmap[ind] ]] = int(line[1])

  alldata['L'] = L
