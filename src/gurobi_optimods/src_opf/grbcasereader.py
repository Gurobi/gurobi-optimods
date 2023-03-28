import math
import cmath
import time

from .myutils import *
from .log import *


class Bus:
    """Bus class"""

    # This is a class that describes a bus in a power system, including loads,
    # voltage limits, type of bus (generator, reference, etc.), ID and count.
    # The data structure also keeps track of branches incident with the bus
    # See comments below

    def __init__(
        self, count, nodeID, nodetype, Pd, Qd, Gs, Bs, Vbase, Vmax, Vmin, busline0
    ):
        self.genidsbycount = []  # array of generator IDs at this bus
        self.frombranchids = {}  # branches where this bus is the 'from' bus
        self.tobranchids = {}  # branches where this bus is the 'to' bus
        self.count = count  # bus count
        self.nodeID = nodeID  # ID of bus
        self.nodetype = nodetype
        self.Pd = Pd  # active load
        self.Qd = Qd  # reactive load
        self.Gs = Gs  # shunt admittance parameter
        self.Bs = Bs  # shunt admittance parameter
        self.Vbase = Vbase  # voltage base
        self.Vmax = Vmax
        self.Vmin = Vmin
        self.busline0 = busline0  # line location of bus within input file
        self.inputvoltage = False
        self.cffvarind = -1  # the following four fields are variable indices
        self.Pinjvarind = -1
        self.Qinjvarind = -1
        self.vvarind = -1
        self.thetavarind = -1
        self.Pbalance = 0
        self.inputV = 0  # input voltage
        self.inpute = 0  # e-value (input or derived from input voltage solution)
        self.inputf = 0  # e-value (input or derived from input voltage solution)
        self.outdegree = self.indegree = self.degree = 0

        self.lat = -1
        self.lon = -1

    def getbusline0(self):
        return self.busline0

    def addgenerator(self, log, generatorcount1, generator):
        self.genidsbycount.append(generatorcount1)
        loud = False
        if loud:
            log.joint(
                " added generator # "
                + str(generatorcount1)
                + " to bus ID "
                + str(self.nodeID)
            )
            log.joint(
                " Pmax " + str(generator.Pmax) + " Pmin " + str(generator.Pmin) + "\n"
            )

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
    """Branch class"""

    # Branch class

    def __init__(
        self,
        log,
        count,
        f,
        count_f,
        t,
        count_t,
        r,
        x,
        bc,
        rateAmva,
        rateBmva,
        rateCmva,
        ratio,
        angle,
        maxangle,
        minangle,
        status,
        defaultlimit,
        branchline0,
    ):
        self.count = count
        self.f = f  # bus ID for from bus
        self.t = t  # bus ID for to bus
        self.count_f = count_f  # count for from bus
        self.count_t = count_t  # count for to bus
        self.r = r  # resistance
        self.x = x  # reactance
        self.bc = bc  # branch charging admittance
        self.count = count
        self.branchline0 = branchline0
        self.rateAmva = rateAmva  # the following three parameters
        self.rateBmva = rateBmva  # describe branch limits
        self.limit = rateAmva
        self.unboundedlimit = False
        self.constrainedflow = 1

        if self.limit == 0:
            self.limit = defaultlimit
            self.constrainedflow = 0
            self.unboundedlimit = True

        if ratio == 0:
            ratio = 1

        self.ratio = ratio
        self.angle = angle
        self.angle_rad = math.pi * angle / 180.0
        self.maxangle = maxangle
        self.maxangle_rad = math.pi * maxangle / 180.0
        self.minangle = minangle
        self.minangle_rad = math.pi * minangle / 180.0
        self.upperanglenone = 0

        if maxangle == 360 or maxangle == 0:
            self.maxangle_rad = 2 * math.pi
            self.upperanglenone = 1

        self.loweranglenone = 0

        if minangle == -360 or minangle == 0:
            self.minangle_rad = -2 * math.pi
            self.loweranglenone = 1

        self.invratio2 = invratio2 = 1 / ratio**2
        self.multtf = multtf = 1 / (ratio * cmath.exp(1j * self.angle_rad))
        self.multft = multft = 1 / (ratio * cmath.exp(-1j * self.angle_rad))
        # print 'multtf', multtf
        self.status = status
        self.z = z = r + x * 1j
        self.y = y = 1 / z
        self.ynorm2 = (y.real) * (y.real) + (y.imag) * (y.imag)
        self.Yff = Yff = (y + bc / 2 * 1j) * invratio2
        self.Yft = Yft = -y * multft
        self.Ytf = Ytf = -y * multtf
        self.Ytt = Ytt = y + bc / 2 * 1j
        self.Gff = Gff = (self.Yff).real
        self.Bff = Bff = (self.Yff).imag
        self.Gft = Gft = (self.Yft).real
        self.Bft = Bft = (self.Yft).imag
        self.Gtf = Gtf = (self.Ytf).real
        self.Btf = Btf = (self.Ytf).imag
        self.Gtt = Gtt = (self.Ytt).real
        self.Btt = Btt = (self.Ytt).imag

        self.isacline = (ratio == 1) and (self.angle_rad == 0)
        self.nongaining = (
            (self.Gff >= 0)
            and (self.Gtt >= 0)
            and (self.Gff >= -self.Gft)
            and (self.Gtt >= -self.Gtf)
        )

        self.inputcs = False
        self.inputc = 2
        self.inputs = 2
        self.cftvarind = -1
        self.sftvarind = -1
        self.Pftvarind = -1
        self.Qftvarind = -1
        self.switchvarind = -1

        loud = False
        if loud:
            log.joint("\nbr " + str(count) + " f " + str(f) + " t " + str(t) + "\n")
            log.joint("   idf " + str(count_f) + " idt " + str(count_t) + "\n")
            log.joint("   r " + str(r) + " x " + str(x) + " bb " + str(bc) + "\n")
            log.joint(
                "   ratio "
                + str(self.ratio)
                + " angle "
                + str(angle)
                + " angle_rad: "
                + str(self.angle_rad)
                + "\n"
            )
            log.joint("   y " + str(y) + "\n")
            log.joint(
                "       Yff "
                + str(Yff)
                + " , Yft "
                + str(Yft)
                + " , Ytf "
                + str(Ytf)
                + " , Ytt "
                + str(Ytt)
                + "\n"
            )

    def getbranchline0(self):
        return self.branchline0

    def show(self, log):
        log.joint(" < " + str(self.f) + " , " + str(self.t) + " > ")
        log.joint(" r " + str(self.r) + " x " + str(self.x) + " bc " + str(self.bc))
        log.joint("\n")
        log.joint(" ra " + str(self.ratio) + " ang " + str(self.angle))
        log.joint("\n")


class Gen:
    """Gen class"""

    # Generator class
    # See Matpower manual for more details

    def __init__(self, count, nodeID, Pg, Qg, status, Pmax, Pmin, Qmax, Qmin, line0):
        self.count = count
        self.nodeID = nodeID  # ID of bus holding gen
        self.Pg = Pg  # active power generation in existing solution
        self.Qg = Qg  # reactive power generation in existing solution
        self.status = status  # on or off
        self.Pmax = Pmax  # max, min active and reactive power limits
        self.Pmin = Pmin
        self.Qmax = Qmax
        self.Qmin = Qmin
        self.line0 = line0
        self.costlinenum = -1
        self.Pvarind = -1
        self.Qvarind = -1

    def addcost(self, log, costvector, linenum):
        self.costvector = costvector
        self.costdegree = len(costvector) - 1
        self.costlinenum = linenum

    def getline0(self):
        return self.line0


def read_case(alldata):
    """Read case file and fill data dictionary"""

    log = alldata["log"]
    casefilename = alldata["casefilename"]
    starttime = time.time()

    log.joint("Reading case file %s\n" % casefilename)
    f = open(casefilename, "r")
    lines = f.readlines()
    f.close()

    # Read all lines of the casefile
    read_case_thrulines(log, alldata, lines)
    alldata["casefilelines"] = lines

    endtime = time.time()

    log.joint("Reading time: %f s\n" % (endtime - starttime))


def read_case_thrulines(log, alldata, lines):
    """Read thru all lines of case file and fill data dictionary"""

    numlines = len(lines)
    lookingforbus = 1
    linenum = 2
    numisolated = 0
    baseMVA = 0.0
    alldata["refbus"] = -1

    # Read file line by line
    while linenum <= numlines:
        line = lines[linenum - 1]
        thisline = line.split()

        # Skip empty line
        if len(thisline) <= 0:
            linenum += 1
            continue

        # Skip unnecessary lines
        theword = thisline[0]
        if len(theword) < 4 or theword[0:4] != "mpc.":
            linenum += 1
            continue

        log.joint("  Found %s on line %d\n" % (theword, linenum))
        if theword == "mpc.baseMVA":
            tmp = thisline[2]
            # Trim ; if present
            if tmp[len(tmp) - 1] == ";":
                tmp = tmp[: len(tmp) - 1]

            alldata["baseMVA"] = float(tmp)
            baseMVA = alldata["baseMVA"]
            log.joint("    baseMVA: %f\n" % baseMVA)

        elif theword == "mpc.bus":
            buses = {}
            IDtoCountmap = {}
            lookingforendofbus = 1
            slackbus = -1
            numbuses = 0
            numPload = 0
            sumPd = sumQd = 0
            linenum += 1
            # Read bus section
            while lookingforendofbus and linenum <= numlines:
                line = lines[linenum - 1]
                thisline = line.split()
                length = len(thisline)

                # Look for end of section
                if thisline[0] == "];":
                    log.joint("    Found end of bus section on line %d\n" % linenum)
                    lookingforendofbus = 0
                    break

                numbuses += 1
                if thisline[1] == "3":
                    slackbus = int(thisline[0])
                    log.joint("    Slack bus: %d\n" % slackbus)

                if thisline[0] != "%":
                    nodeID = int(thisline[0])
                    nodetype = int(thisline[1])

                    if (
                        nodetype != 1
                        and nodetype != 2
                        and nodetype != 3
                        and nodetype != 4
                    ):
                        raise ValueError(
                            "Bad bus %s has type %s." % (thisline[0], thisline[1])
                        )

                    if nodetype == 4:
                        numisolated += 1

                    # Trim ; if present
                    Vmin = thisline[12]
                    if Vmin[len(Vmin) - 1] == ";":
                        Vmin = Vmin[: len(Vmin) - 1]

                    Vmin = float(Vmin)
                    Pd = float(thisline[2])
                    Qd = float(thisline[3])
                    Gs = float(thisline[4])
                    Bs = float(thisline[5])
                    Vbase = float(thisline[9])
                    Vmax = float(thisline[11])

                    if baseMVA == 0.0:
                        raise ValueError("baseMVA not available before bus section.")

                    buses[numbuses] = Bus(
                        numbuses,
                        nodeID,
                        nodetype,
                        Pd / baseMVA,
                        Qd / baseMVA,
                        Gs / baseMVA,
                        Bs / baseMVA,
                        Vbase,
                        Vmax,
                        Vmin,
                        linenum - 1,
                    )

                    if nodetype == 3:
                        log.joint(
                            "    Bus %d ID %d is the reference bus\n"
                            % (numbuses, nodeID)
                        )
                        alldata["refbus"] = numbuses

                    if nodetype == 1 or nodetype == 2 or nodetype == 3:
                        sumPd += Pd
                        sumQd += Qd

                    IDtoCountmap[nodeID] = numbuses
                    numPload += Pd > 0
                else:
                    nodeID = int(thisline[1])
                    nodetype = int(thisline[2])
                    log.joint(
                        "bus %d nodeID %d is isolated and has type %d\n"
                        % (numbuses, nodeID, nodetype)
                    )
                    log.joint("   setting it to type 4\n")
                    nodetype = 4
                    break_exit("isolated")

                    if (
                        nodetype != 1
                        and nodetype != 2
                        and nodetype != 3
                        and nodetype != 4
                    ):
                        raise ValueError(
                            "Bad bus %s has type %s." % (thisline[0], thisline[1])
                        )

                    if nodetype == 4:
                        numisolated += 1

                    # Trim ; if present
                    Vmin = 0

                    Pd = 0
                    Qd = 0
                    Gs = 0
                    Bs = 0
                    Vbase = float(thisline[10])
                    Vmax = 0

                    if baseMVA == 0.0:
                        raise ValueError("baseMVA not available before bus section.")

                    buses[numbuses] = Bus(
                        numbuses,
                        nodeID,
                        nodetype,
                        Pd / baseMVA,
                        Qd / baseMVA,
                        Gs / baseMVA,
                        Bs / baseMVA,
                        Vbase,
                        Vmax,
                        Vmin,
                        linenum - 1,
                    )

                    if nodetype == 3:
                        log.joint(
                            "    Bus %d ID %d is the reference bus\n"
                            % (numbuses, nodeID)
                        )
                        alldata["refbus"] = numbuses

                    if nodetype == 1 or nodetype == 2 or nodetype == 3:
                        sumPd += Pd
                        sumQd += Qd

                    IDtoCountmap[nodeID] = numbuses
                    numPload += Pd > 0

                linenum += 1
            # Finished reading bus section
            alldata["buses"] = buses
            alldata["numbuses"] = numbuses
            alldata["sumPd"] = sumPd
            alldata["sumQd"] = sumQd
            alldata["IDtoCountmap"] = IDtoCountmap
            alldata["slackbus"] = slackbus

            log.joint("    sumloadPd %f numPload %d\n" % (sumPd, numPload))
            log.joint("    sumloadQd %f\n" % sumQd)

            if lookingforendofbus:
                raise ValueError("Could not find bus data section.")

            if slackbus < 0:
                log.joint("    Could not find slack bus\n")

            log.joint("    %d buses\n" % numbuses)
            if numisolated > 0:
                log.joint("    isolated buses: %d\n" % numisolated)

        elif theword == "mpc.gen":
            gens = {}
            lookingforendofgen = 1
            gencount1 = 0
            summaxgenP = summaxgenQ = 0
            linenum += 1

            # Read gen section
            while lookingforendofgen and linenum <= numlines:
                line = lines[linenum - 1]
                thisline = line.split()

                # Look for end of section
                if thisline[0] == "];":
                    alldata["endofgen"] = linenum
                    log.joint("    Found end of gen section on line %d\n" % linenum)
                    lookingforendofgen = 0
                    break

                gencount1 += 1

                nodeID = int(thisline[0])
                Pg = float(thisline[1])
                Qg = float(thisline[2])
                status = int(thisline[7])
                Pmax = float(thisline[8])
                Pmin = float(thisline[9])
                Qmax = float(thisline[3])
                Qmin = float(thisline[4])

                if status <= 0:
                    status = 0
                else:
                    status = 1

                if baseMVA == 0.0:
                    raise ValueError("baseMVA not available before gen section.")

                if nodeID in IDtoCountmap.keys():
                    idgencount1 = IDtoCountmap[nodeID]
                    gens[gencount1] = Gen(
                        gencount1,
                        nodeID,
                        Pg,
                        Qg,
                        status,
                        Pmax / baseMVA,
                        Pmin / baseMVA,
                        Qmax / baseMVA,
                        Qmin / baseMVA,
                        linenum - 1,
                    )
                    buses[idgencount1].addgenerator(log, gencount1, gens[gencount1])

                    if (
                        buses[idgencount1].nodetype == 2
                        or buses[idgencount1].nodetype == 3
                    ):  # But not 4
                        summaxgenP += Pmax
                        summaxgenQ += Qmax

                else:
                    raise ValueError(
                        "Generator # %d in nonexistent bus ID %d." % (gencount1, nodeID)
                    )

                linenum += 1
            # Finished reading gen section
            if lookingforendofgen:
                raise ValueError("Could not find end of generator section.")

            alldata["gens"] = gens
            alldata["numgens"] = len(gens)
            busgencount = 0

            for bus in buses.values():
                busgencount += len(bus.genidsbycount) > 0

            alldata["busgencount"] = busgencount
            alldata["summaxgenP"] = summaxgenP
            alldata["summaxgenQ"] = summaxgenQ

            log.joint("    number of generators: %d\n" % alldata["numgens"])
            log.joint("    number of buses with gens: %d\n" % busgencount)
            log.joint("    summaxPg %f summaxQg %f\n" % (summaxgenP, summaxgenQ))

        elif theword == "mpc.branch":
            branches = {}
            defaultlimit = 1e20
            lookingforendofbranch = 1
            numbranches = 0
            activebranches = 0
            zerolimit = 0
            linenum += 1

            # Read branch section
            while lookingforendofbranch and linenum <= numlines:
                line = lines[linenum - 1]
                thisline = line.split()

                # Look for end of branch section
                if thisline[0] == "];":
                    log.joint("    Found end of branch section on line %d\n" % linenum)
                    lookingforendofbranch = 0
                    break

                numbranches += 1
                f = int(thisline[0])
                t = int(thisline[1])
                r = float(thisline[2])
                x = float(thisline[3])
                bc = float(thisline[4])
                rateA = float(thisline[5])
                rateB = float(thisline[6])
                rateC = float(thisline[7])
                ratio = float(thisline[8])
                if ratio == 0:
                    ratio = 1.0
                angle = float(thisline[9])
                status = int(thisline[10])
                minangle = float(thisline[11])
                maxangle = thisline[12]
                # Trim ; at the end of line
                if maxangle[len(maxangle) - 1] == ";":
                    maxangle = maxangle[: len(maxangle) - 1]

                maxangle = float(maxangle)

                if maxangle < minangle:
                    raise ValueError(
                        "Branch # %d has illegal angle constraints\n" % numbranches
                    )

                count_f = IDtoCountmap[f]
                count_t = IDtoCountmap[t]

                if baseMVA == 0.0:
                    raise ValueError("baseMVA not available before branch section.")

                if False:  # rateA < 1e-10: # TODO do we need it?
                    log.joint(
                        "Warning. Branch %d from %d to %d has unbounded rateA.\n"
                        % (numbranches, f, t)
                    )

                if True:  # TODO do we need the if?
                    branches[numbranches] = Branch(
                        log,
                        numbranches,
                        f,
                        count_f,
                        t,
                        count_t,
                        r,
                        x,
                        bc,
                        rateA / baseMVA,
                        rateB / baseMVA,
                        rateC / baseMVA,
                        ratio,
                        angle,
                        maxangle,
                        minangle,
                        status,
                        defaultlimit,
                        linenum - 1,
                    )
                    zerolimit += branches[numbranches].constrainedflow == 0
                    activebranches += 1
                    buses[count_f].addfrombranch(log, numbranches)
                    buses[count_t].addtobranch(log, numbranches)

                linenum += 1

            # Finished reading branch section
            if lookingforendofbranch:
                raise ValueError("Could not find end of branch section.")

            alldata["branches"] = branches
            alldata["numbranches"] = numbranches
            log.joint(
                "    numbranches: %d active: %d\n" % (numbranches, activebranches)
            )
            if zerolimit > 0:
                log.joint("    ---> %d unconstrained.\n" % zerolimit)

        elif theword == "mpc.gencost":
            lookingforendofgencost = 1
            gencostcount = 1
            linenum += 1

            # Read gen cost section
            while lookingforendofgencost and linenum <= numlines:
                line = lines[linenum - 1]
                thisline = line.split()

                # Look for end of gen section
                if thisline[0] == "];":
                    log.joint("    Found end of gencost section on line %d\n" % linenum)
                    alldata["endofgencost"] = linenum
                    lookingforendofgencost = 0
                    break

                if gencostcount <= gencount1:
                    costtype = int(thisline[0])
                    if costtype != 2:
                        raise ValueError(
                            "Cost of generator %d is not polynomial\n" % gencostcount
                        )

                    degree = int(thisline[3]) - 1
                    skip1 = 0

                    if degree == 3:
                        coeff = float(thisline[4])
                        skip1 = 1
                        if coeff != 0:
                            raise ValueError(
                                "Degree of cost function for generator %d equals %d which is illegal\n"
                                % (gencostcount, degree)
                            )
                        else:
                            degree = 2

                    elif degree > 2 or degree < 0:
                        raise ValueError(
                            "Degree of cost function for generator %d equals %d which is illegal\n"
                            % (gencostcount, degree)
                        )

                    costvector = [0 for j in range(degree + 1)]

                    for j in range(degree + 1):
                        tmp = thisline[4 + skip1 + j]
                        # print(j, tmp)

                        if tmp[len(tmp) - 1] == ";":
                            tmp = tmp[: len(tmp) - 1]

                        costvector[j] = float(tmp)
                        costvector[j] *= (baseMVA) ** (degree - j)

                    gens[gencostcount].addcost(log, costvector, linenum)

                    # break_exit('gen')

                else:
                    raise ValueError(
                        "Read %d gen costs but only %d generators\n"
                        % (gencostcount, gencount1)
                    )

                gencostcount += 1
                linenum += 1
            # Finished reading gen cost section
            if lookingforendofgencost:
                raise ValueError("Could not find end of gencost section.")

            linenum += 1

        linenum += 1
    # Finished reading file line by line


def readvoltsfile(log, alldata):
    """Read volts file and fill data dictionary"""

    voltsfilename = alldata["voltsfilename"]

    log.joint("  Reading volts file %s\n" % voltsfilename)
    f = open(voltsfilename, "r")
    lines = f.readlines()
    f.close()

    inputvolts = {}
    numread = 0

    for linenum in range(len(lines)):
        thisline = lines[linenum].split()

        if len(thisline) <= 0:
            continue

        if thisline[0] == "bus":
            angle_rad = float(thisline[5]) * math.pi / 180
            busid = int(thisline[1])
            inputvolts[busid] = (float(thisline[3]), angle_rad)
            numread += 1

        elif thisline[0] == "END":
            break

        else:
            raise ValueError(
                "Illegal input %s on line %s in volts file\n"
                % (thisline[0], lines[linenum])
            )

    log.joint("Read %d input voltages\n" % numread)
    alldata["inputvolts"] = inputvolts


def readflowsfile(log, alldata):
    """Read flows file and fill data dictionary"""

    flowsfilename = alldata["flowsfilename"]

    log.joint("  Reading flows file %s\n" % flowsfilename)
    f = open(flowsfilename, "r")
    lines = f.readlines()
    f.close()

    baseMVA = alldata["baseMVA"]
    inputPf = {}
    inputQf = {}
    inputPt = {}
    inputQt = {}
    numread = 0

    for linenum in range(len(lines)):
        thisline = lines[linenum].split()

        if len(thisline) <= 0:
            continue

        if thisline[0] == "branch":
            branchid = int(thisline[1])
            inputPf[branchid] = float(thisline[7]) / baseMVA
            inputPt[branchid] = float(thisline[9]) / baseMVA
            inputQf[branchid] = float(thisline[11]) / baseMVA
            inputQt[branchid] = float(thisline[13]) / baseMVA

            numread += 1
        elif thisline[0] == "END":
            break

        else:
            raise ValueError(
                "Illegal input %s on line %s in flows file\n"
                % (thisline[0], lines[linenum])
            )

    log.joint("Read %d input flows\n" % numread)

    alldata["inputPf"] = inputPf
    alldata["inputPt"] = inputPt
    alldata["inputQf"] = inputQf
    alldata["inputQt"] = inputQt


def writegv(log, alldata, gvfilename):
    """Write gv file needed for graphical image"""

    f = open(gvfilename, "w")
    log.joint("Writing to gv file %s\n" % gvfilename)

    f.write("graph {\n")

    for bus in alldata["buses"].values():
        f.write("     " + str(bus.nodeID) + ";\n")

    for branch in alldata["branches"].values():
        f.write("     " + str(branch.f) + " -- " + str(branch.t) + ";\n")

    f.write("}\n")
    f.close()


def generateinputcs(log, alldata):
    """Description"""
    # Generates values for the Jabr c and s variables from input voltages

    log.joint("  Generating input c,s values\n")

    inputcc = {}
    inputcs = {}
    inputvolts = alldata["inputvolts"]
    branches = alldata["branches"]
    IDtoCountmap = alldata["IDtoCountmap"]

    for busid in inputvolts:
        M = inputvolts[busid][0]
        inputcc[busid] = M * M

    for branch in branches.values():
        f = branch.f
        t = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]

        if f in inputvolts and t in inputvolts:
            Mf = inputvolts[f][0]
            Mt = inputvolts[t][0]
            af = inputvolts[f][1]
            at = inputvolts[t][1]
            angle = af - at
            inputcs[f, t] = (Mf * Mt * math.cos(angle), Mf * Mt * math.sin(angle))
            inputcs[t, f] = (Mf * Mt * math.cos(angle), -Mf * Mt * math.sin(angle))

    alldata["inputcc"] = inputcc
    alldata["inputcs"] = inputcs


def generateinputeandf(log, alldata):
    """Description"""
    # Generates rectangular coordinates for voltages from input solution (given in polar coordinates)

    log.joint("  generating input e,f values\n")

    inputve = {}
    inputvf = {}
    inputvolts = alldata["inputvolts"]

    for busid in inputvolts:
        M = inputvolts[busid][0]
        A = inputvolts[busid][1]
        inputve[busid] = M * math.cos(A)
        inputvf[busid] = M * math.sin(A)

    alldata["inputve"] = inputve
    alldata["inputvf"] = inputvf
