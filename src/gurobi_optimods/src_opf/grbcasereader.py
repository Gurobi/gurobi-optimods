import math
import cmath
import time
import logging


class Bus:
    """
    Bus class
    This is a class that describes a bus in a power system, including loads,
    voltage limits, type of bus (generator, reference, etc.), ID and count.
    The data structure also keeps track of branches incident with the bus
    """

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

    def addgenerator(self, generatorcount1, generator):
        logger = logging.getLogger("OpfLogger")
        self.genidsbycount.append(generatorcount1)
        loud = False
        if loud:
            logger.info(
                " added generator # "
                + str(generatorcount1)
                + " to bus ID "
                + str(self.nodeID)
            )
            logger.info(" Pmax " + str(generator.Pmax) + " Pmin " + str(generator.Pmin))

    def addfrombranch(self, id):
        quant = len(self.frombranchids)
        self.frombranchids[quant] = id
        self.outdegree += 1
        self.degree += 1

    def addtobranch(self, id):
        quant = len(self.tobranchids)
        self.tobranchids[quant] = id
        self.indegree += 1
        self.degree += 1


class Branch:
    """
    Branch class
    TODO-Dan Please add more description if needed
    """

    def __init__(
        self,
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

        loud = False  # TODO-Dan Do we need this?
        if loud:
            logger = logging.getLogger("OpfLogger")
            logger.info("\nbr " + str(count) + " f " + str(f) + " t " + str(t))
            logger.info("   idf " + str(count_f) + " idt " + str(count_t))
            logger.info("   r " + str(r) + " x " + str(x) + " bb " + str(bc))
            logger.info(
                "   ratio "
                + str(self.ratio)
                + " angle "
                + str(angle)
                + " angle_rad: "
                + str(self.angle_rad)
            )
            logger.info("   y " + str(y))
            logger.info(
                "       Yff "
                + str(Yff)
                + " , Yft "
                + str(Yft)
                + " , Ytf "
                + str(Ytf)
                + " , Ytt "
                + str(Ytt)
            )

    def getbranchline0(self):
        return self.branchline0

    def show(self):
        logger = logging.getLogger("OpfLogger")
        logger.info(" < " + str(self.f) + " , " + str(self.t) + " > ")
        logger.info(" r " + str(self.r) + " x " + str(self.x) + " bc " + str(self.bc))
        logger.info(" ra " + str(self.ratio) + " ang " + str(self.angle))


class Gen:
    """
    Generator class
    See Matpower manual for more details
    TODO-Dan Could you please add a link
    """

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

    def addcost(self, costvector, linenum):
        self.costvector = costvector
        self.costdegree = len(costvector) - 1
        self.costlinenum = linenum

    def addcost_plain(self, costvector):
        self.costvector = costvector
        self.costdegree = len(costvector) - 1
        self.costlinenum = -1

    def showcostvector(self):
        logger = logging.getLogger("OpfLogger")
        logger.info(self.costvector)
        for i in range(0, self.costdegree + 1):
            logger.info(i, self.costvector[i], " ", end="")

    def getline0(self):
        return self.line0


def read_case(alldata, case):
    if type(case) is dict:
        # Fill alldata dictionary from case dictionary
        build_data_struct(alldata, case)
    else:
        # casefilename has been set in read_settings
        if alldata["casefilename"] == None:
            raise ValueError("No casefile provided.")

        # Read case file holding OPF network data.
        read_casefile(alldata)


def build_data_struct(alldata, case_dict):
    """
    Fills alldata dictionary out of a given casefile dictionary
    """

    logger = logging.getLogger("OpfLogger")
    buses = {}
    IDtoCountmap = {}
    slackbus = -1
    numbuses = 0
    numPload = 0
    sumPd = sumQd = 0
    numisolated = 0

    logger.info("Building main data structure from dictionary.")

    dict_buses = case_dict["buses"]
    baseMVA = alldata["baseMVA"] = case_dict["baseMVA"]

    logger.info("Buses.")

    for dbus in dict_buses.values():
        numbuses += 1
        count = dbus["count"]
        nodeID = dbus["nodeID"]
        nodetype = dbus["nodetype"]

        if nodetype != 1 and nodetype != 2 and nodetype != 3 and nodetype != 4:
            raise ValueError("Bad bus %s has type %s." % (count, nodetype))

        if nodetype == 3:
            slackbus = count
            logger.info("    Bus %d ID %d is the reference bus." % (count, nodeID))
            alldata["refbus"] = count

        if nodetype == 4:
            numisolated += 1

        Vmin = dbus["Vmin"]
        Pd = dbus["Pd"]
        Qd = dbus["Qd"]
        Gs = dbus["Gs"]
        Bs = dbus["Bs"]
        Vbase = dbus["Vbase"]
        Vmax = dbus["Vmax"]
        Vmin = dbus["Vmin"]
        lnum = dbus["lnum"]

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
            lnum,
        )

        if nodetype == 1 or nodetype == 2 or nodetype == 3:
            sumPd += Pd
            sumQd += Qd

        IDtoCountmap[nodeID] = count
        numPload += Pd > 0

    alldata["buses"] = buses
    alldata["numbuses"] = numbuses
    alldata["sumPd"] = sumPd
    alldata["sumQd"] = sumQd
    alldata["IDtoCountmap"] = IDtoCountmap
    alldata["slackbus"] = slackbus

    logger.info("    sumloadPd %f numPload %d" % (sumPd, numPload))
    logger.info("    sumloadQd %f" % sumQd)

    if slackbus < 0:
        logger.info("    Could not find slack bus.")

    logger.info("    %d buses" % numbuses)
    if numisolated > 0:
        logger.info("    isolated buses: %d" % numisolated)

    logger.info("Branches.")
    branches = {}
    defaultlimit = 1e20
    numbranches = 0
    activebranches = 0
    zerolimit = 0

    dict_branches = case_dict["branches"]
    branches = {}

    for dbranch in dict_branches.values():
        # print(dbranch)
        numbranches += 1
        brcnt1 = dbranch["branchcount1"]
        f = dbranch["f"]
        t = dbranch["t"]
        r = dbranch["r"]
        x = dbranch["x"]
        bc = dbranch["bc"]
        rateA = dbranch["rateA"]
        rateB = dbranch["rateB"]
        rateC = dbranch["rateC"]
        ratio = dbranch["ratio"]
        if ratio == 0:
            ratio = 1.0  # to be sure
        angle = dbranch["angle"]
        status = dbranch["status"]
        minangle = dbranch["minangle"]
        maxangle = dbranch["maxangle"]
        if maxangle < minangle:
            raise ValueError("Branch # %d has illegal angle constraints." % numbranches)

        lnum = dbranch["lnum"]

        count_f = IDtoCountmap[f]
        count_t = IDtoCountmap[t]
        branches[numbranches] = Branch(
            brcnt1,
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
            lnum,
        )
        activebranches += 1
        buses[count_f].addfrombranch(brcnt1)
        buses[count_t].addtobranch(brcnt1)

    alldata["branches"] = branches
    alldata["numbranches"] = numbranches
    logger.info("    numbranches: %d active: %d" % (numbranches, activebranches))
    if zerolimit > 0:
        logger.info("    ---> %d unconstrained" % zerolimit)

    logger.info("Generators.")
    gens = {}
    dict_gens = case_dict["generators"]
    summaxgenP = summaxgenQ = 0

    for dgen in dict_gens.values():
        gencount1 = dgen["gencount1"]
        nodeID = dgen["nodeID"]
        Pg = dgen["Pg"]
        Qg = dgen["Qg"]
        status = dgen["status"]
        Pmax = dgen["Pmax"]
        Pmin = dgen["Pmin"]
        Qmax = dgen["Qmax"]
        Qmin = dgen["Qmin"]
        lnum = dgen["lnum"]

        idgencount1 = -1
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
                lnum,
            )

            buses[idgencount1].addgenerator(gencount1, gens[gencount1])

            if (
                buses[idgencount1].nodetype == 2 or buses[idgencount1].nodetype == 3
            ):  # But not 4
                summaxgenP += Pmax
                summaxgenQ += Qmax
        else:
            raise ValueError(
                "Generator # %d in nonexistent bus ID %d." % (gencount1, nodeID)
            )

    alldata["gens"] = gens
    alldata["numgens"] = len(gens)
    busgencount = 0

    for bus in buses.values():
        busgencount += len(bus.genidsbycount) > 0

    alldata["busgencount"] = busgencount
    alldata["summaxgenP"] = summaxgenP
    alldata["summaxgenQ"] = summaxgenQ

    logger.info("    number of generators: %d" % alldata["numgens"])
    logger.info("    number of buses with gens: %d" % busgencount)
    logger.info("    summaxPg %f summaxQg %f" % (summaxgenP, summaxgenQ))

    logger.info("Generator cost vectors.")
    gencoststruct = case_dict["generator_cost_structure"]

    for count in range(1, alldata["numgens"] + 1):
        gens[count].addcost_plain(gencoststruct[count]["costvector"])


def read_case_build_dict(alldata):
    """
    Read case file and construct a dict out of it.

    This is a helper function for translating casefiles into dict format.
    """

    logger = logging.getLogger("OpfLogger")
    casefilename = alldata["casefilename"]
    starttime = time.time()

    logger.info("Reading case file %s." % casefilename)
    f = open(casefilename, "r")
    lines = f.readlines()
    f.close()

    logger.info("Now building case dictionary.")

    # Read all lines of the casefile
    alldata["casefilelines"] = lines

    case_dict = read_case_build_dict_thrulines(lines)

    endtime = time.time()

    logger.info("Reading time: %f s." % (endtime - starttime))

    build_data_struct(alldata, case_dict)

    return case_dict


def read_case_build_dict_thrulines(lines):
    """
    Read thru all lines of case file and fill data dictionary.
    This is a helper function for translating casefiles into dict format.
    """

    logger = logging.getLogger("OpfLogger")
    numlines = len(lines)
    lookingforbus = 1
    linenum = 2
    baseMVA = 0.0

    case_dict = {}
    case_dict["refbus"] = -1

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

        logger.info("  Found %s on line %d." % (theword, linenum))
        if theword == "mpc.baseMVA":
            tmp = thisline[2]
            # Trim ; if present
            if tmp[len(tmp) - 1] == ";":
                tmp = tmp[: len(tmp) - 1]

            case_dict["baseMVA"] = baseMVA = float(tmp)
            logger.info("    baseMVA: %f" % baseMVA)

        elif theword == "mpc.bus":
            buses = {}
            case_dict["buses"] = buses
            lookingforendofbus = 1
            slackbus = -1
            numbuses = 0
            linenum += 1
            # Read bus section
            while lookingforendofbus and linenum <= numlines:
                line = lines[linenum - 1]
                thisline = line.split()
                length = len(thisline)

                # Look for end of section
                if thisline[0] == "];":
                    logger.info("    Found end of bus section on line %d." % linenum)
                    lookingforendofbus = 0
                    break

                numbuses += 1
                if thisline[1] == "3":
                    slackbus = int(thisline[0])
                    logger.info("    Slack bus: %d" % slackbus)

                if thisline[0] != "%":
                    nodeID = int(thisline[0])
                    nodetype = int(thisline[1])

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

                    buses[numbuses] = {}
                    buses[numbuses]["count"] = numbuses
                    buses[numbuses]["nodeID"] = nodeID
                    buses[numbuses]["nodetype"] = nodetype
                    buses[numbuses]["Pd"] = Pd
                    buses[numbuses]["Qd"] = Qd
                    buses[numbuses]["Gs"] = Gs
                    buses[numbuses]["Bs"] = Bs
                    buses[numbuses]["Vbase"] = Vbase
                    buses[numbuses]["Vmax"] = Vmax
                    buses[numbuses]["Vmin"] = Vmin
                    buses[numbuses]["lnum"] = linenum - 1

                    if nodetype == 3:
                        logger.info(
                            "    Bus %d ID %d is the reference bus."
                            % (numbuses, nodeID)
                        )
                        case_dict["refbus"] = numbuses

                else:
                    nodeID = int(thisline[1])
                    nodetype = int(thisline[2])
                    logger.info(
                        "bus %d nodeID %d is isolated and has type %d."
                        % (numbuses, nodeID, nodetype)
                    )
                    logger.info("   setting it to type 4.")
                    nodetype = 4

                    # Trim ; if present
                    Vmin = 0
                    Pd = 0
                    Qd = 0
                    Gs = 0
                    Bs = 0
                    Vbase = float(thisline[10])
                    Vmax = Vmin = 0

                    buses[numbuses] = {}
                    buses[numbuses]["count"] = numbuses
                    buses[numbuses]["nodeID"] = nodeID
                    buses[numbuses]["nodetype"] = nodetype
                    buses[numbuses]["Pd"] = Pd
                    buses[numbuses]["Qd"] = Qd
                    buses[numbuses]["Gs"] = Gs
                    buses[numbuses]["Bs"] = Bs
                    buses[numbuses]["Vbase"] = Vbase
                    buses[numbuses]["Vmax"] = Vmax
                    buses[numbuses]["Vmin"] = Vmin

                linenum += 1
            # Finished reading bus section
            case_dict["slackbus"] = slackbus

            if lookingforendofbus:
                raise ValueError("Could not find bus data section.")

        elif theword == "mpc.gen":
            case_dict["generators"] = gens = {}
            lookingforendofgen = 1
            gencount1 = 0
            linenum += 1

            # Read gen section
            while lookingforendofgen and linenum <= numlines:
                line = lines[linenum - 1]
                thisline = line.split()

                # Look for end of section
                if thisline[0] == "];":
                    logger.info("    Found end of gen section on line %d." % linenum)
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

                gens[gencount1] = {}
                gens[gencount1]["gencount1"] = gencount1
                gens[gencount1]["nodeID"] = nodeID
                gens[gencount1]["Pg"] = Pg
                gens[gencount1]["Qg"] = Qg
                gens[gencount1]["status"] = status
                gens[gencount1]["Pmax"] = Pmax
                gens[gencount1]["Pmin"] = Pmin
                gens[gencount1]["Qmax"] = Qmax
                gens[gencount1]["Qmin"] = Qmin
                gens[gencount1]["lnum"] = linenum - 1

                linenum += 1
            # Finished reading gen section
            if lookingforendofgen:
                raise ValueError("Could not find end of generator section.")

        elif theword == "mpc.branch":
            case_dict["branches"] = branches = {}
            lookingforendofbranch = 1
            numbranches = 0
            linenum += 1

            # Read branch section
            while lookingforendofbranch and linenum <= numlines:
                line = lines[linenum - 1]
                thisline = line.split()

                # Look for end of branch section
                if thisline[0] == "];":
                    logger.info("    Found end of branch section on line %d." % linenum)
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

                branches[numbranches] = {}
                branches[numbranches]["branchcount1"] = numbranches
                branches[numbranches]["f"] = f
                branches[numbranches]["t"] = t
                branches[numbranches]["r"] = r
                branches[numbranches]["x"] = x
                branches[numbranches]["bc"] = bc
                branches[numbranches]["rateA"] = rateA
                branches[numbranches]["rateB"] = rateB
                branches[numbranches]["rateC"] = rateC
                branches[numbranches]["ratio"] = ratio
                branches[numbranches]["angle"] = angle
                branches[numbranches]["status"] = status
                branches[numbranches]["minangle"] = minangle
                branches[numbranches]["maxangle"] = maxangle
                branches[numbranches]["lnum"] = linenum - 1

                linenum += 1

            # Finished reading branch section
            if lookingforendofbranch:
                raise ValueError("Could not find end of branch section.")

        elif theword == "mpc.gencost":
            case_dict["generator_cost_structure"] = gencoststruct = {}

            lookingforendofgencost = 1
            gencostcount = 1
            linenum += 1

            # Read gen cost section
            while lookingforendofgencost and linenum <= numlines:
                line = lines[linenum - 1]
                thisline = line.split()

                # Look for end of gen section
                if thisline[0] == "];":
                    logger.info(
                        "    Found end of gencost section on line %d." % linenum
                    )
                    lookingforendofgencost = 0
                    break

                if gencostcount <= gencount1:
                    gencoststruct[gencostcount] = {}
                    costtype = int(thisline[0])
                    degree = int(thisline[3]) - 1

                    gencoststruct[gencostcount]["costtype"] = costtype
                    gencoststruct[gencostcount]["degree"] = degree

                    costvector = [0 for j in range(degree + 1)]

                    for j in range(degree + 1):
                        tmp = thisline[4 + j]
                        # print(j, tmp)

                        if tmp[len(tmp) - 1] == ";":
                            tmp = tmp[: len(tmp) - 1]

                        costvector[j] = float(tmp)
                        costvector[j] *= (baseMVA) ** (degree - j)

                    gencoststruct[gencostcount]["costvector"] = costvector

                else:
                    raise ValueError(
                        "Read %d gen costs but only %d generators."
                        % (gencostcount, gencount1)
                    )

                gencostcount += 1
                linenum += 1
            # Finished reading gen cost section
            if lookingforendofgencost:
                raise ValueError("Could not find end of gencost section.")

            case_dict["generator_cost_count"] = gencostcount - 1
            linenum += 1

        linenum += 1
    # Finished reading file line by line

    return case_dict


def read_casefile(alldata):
    # TODO Parse matlab file to matlab data file (MAT) and read it in through Python numpy
    """Read case file and fill data dictionary"""

    logger = logging.getLogger("OpfLogger")
    if alldata["casefilename"] == None:
        raise ValueError("No casefile provided.")

    casefilename = alldata["casefilename"]
    starttime = time.time()

    logger.info("Reading case file %s." % casefilename)
    f = open(casefilename, "r")
    lines = f.readlines()
    f.close()

    # Read all lines of the casefile
    read_case_thrulines(alldata, lines)
    alldata["casefilelines"] = lines

    endtime = time.time()

    logger.info("Reading time: %f s." % (endtime - starttime))


def read_case_thrulines(alldata, lines):
    """Read thru all lines of case file and fill data dictionary"""

    logger = logging.getLogger("OpfLogger")
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

        logger.info("  Found %s on line %d." % (theword, linenum))
        if theword == "mpc.baseMVA":
            tmp = thisline[2]
            # Trim ; if present
            if tmp[len(tmp) - 1] == ";":
                tmp = tmp[: len(tmp) - 1]

            alldata["baseMVA"] = float(tmp)
            baseMVA = alldata["baseMVA"]
            logger.info("    baseMVA: %f." % baseMVA)

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
                    logger.info("    Found end of bus section on line %d." % linenum)
                    lookingforendofbus = 0
                    break

                numbuses += 1
                if thisline[1] == "3":
                    slackbus = int(thisline[0])
                    logger.info("    Slack bus: %d." % slackbus)

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
                        logger.info(
                            "    Bus %d ID %d is the reference bus."
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
                    logger.info(
                        "bus %d nodeID %d is isolated and has type %d."
                        % (numbuses, nodeID, nodetype)
                    )
                    logger.info("   setting it to type 4.")
                    nodetype = 4

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
                        logger.info(
                            "    Bus %d ID %d is the reference bus."
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

            logger.info("    sumloadPd %f numPload %d." % (sumPd, numPload))
            logger.info("    sumloadQd %f." % sumQd)

            if lookingforendofbus:
                raise ValueError("Could not find bus data section.")

            if slackbus < 0:
                logger.info("    Could not find slack bus.")

            logger.info("    %d buses." % numbuses)
            if numisolated > 0:
                logger.info("    isolated buses: %d." % numisolated)

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
                    logger.info("    Found end of gen section on line %d." % linenum)
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
                    buses[idgencount1].addgenerator(gencount1, gens[gencount1])

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

            logger.info("    number of generators: %d." % alldata["numgens"])
            logger.info("    number of buses with gens: %d." % busgencount)
            logger.info("    summaxPg %f summaxQg %f." % (summaxgenP, summaxgenQ))

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
                    logger.info("    Found end of branch section on line %d." % linenum)
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
                        "Branch # %d has illegal angle constraints." % numbranches
                    )

                count_f = IDtoCountmap[f]
                count_t = IDtoCountmap[t]

                if baseMVA == 0.0:
                    raise ValueError("baseMVA not available before branch section.")

                if False:  # rateA < 1e-10: # TODO do we need it?
                    logger.info(
                        "Warning. Branch %d from %d to %d has unbounded rateA."
                        % (numbranches, f, t)
                    )

                if True:  # TODO do we need the if?
                    branches[numbranches] = Branch(
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
                    buses[count_f].addfrombranch(numbranches)
                    buses[count_t].addtobranch(numbranches)

                linenum += 1

            # Finished reading branch section
            if lookingforendofbranch:
                raise ValueError("Could not find end of branch section.")

            alldata["branches"] = branches
            alldata["numbranches"] = numbranches
            logger.info(
                "    numbranches: %d active: %d." % (numbranches, activebranches)
            )
            if zerolimit > 0:
                logger.info("    ---> %d unconstrained." % zerolimit)

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
                    logger.info(
                        "    Found end of gencost section on line %d." % linenum
                    )
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

                    gens[gencostcount].addcost(costvector, linenum)

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


def readvoltsfile(alldata):
    """Read volts file and fill data dictionary"""

    logger = logging.getLogger("OpfLogger")
    if alldata["voltsfilename"] == None:
        raise ValueError("No voltsfilename provided.")

    voltsfilename = alldata["voltsfilename"]

    logger.info("  Reading volts file %s." % voltsfilename)
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

    logger.info("Read %d input voltages." % numread)
    alldata["inputvolts"] = inputvolts


def readflowsfile(alldata):
    """Read flows file and fill data dictionary"""

    logger = logging.getLogger("OpfLogger")
    flowsfilename = alldata["flowsfilename"]

    logger.info("  Reading flows file %s." % flowsfilename)
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
                "Illegal input %s on line %s in flows file."
                % (thisline[0], lines[linenum])
            )

    logger.info("Read %d input flows." % numread)

    alldata["inputPf"] = inputPf
    alldata["inputPt"] = inputPt
    alldata["inputQf"] = inputQf
    alldata["inputQt"] = inputQt


def writegv(alldata, gvfilename):
    """Write gv file needed for graphical image"""

    logger = logging.getLogger("OpfLogger")
    f = open(gvfilename, "w")
    logger.info("Writing to gv file %s." % gvfilename)

    f.write("graph {\n")

    for bus in alldata["buses"].values():
        f.write("     " + str(bus.nodeID) + ";\n")

    for branch in alldata["branches"].values():
        f.write("     " + str(branch.f) + " -- " + str(branch.t) + ";\n")

    f.write("}\n")
    f.close()


def generateinputcs(alldata):
    """
    Generates values for the Jabr c and s variables from input voltages
    """

    logger = logging.getLogger("OpfLogger")
    logger.info("  Generating input c,s values.")

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


def generateinputeandf(alldata):
    """
    Generates rectangular coordinates for voltages from input solution (given in polar coordinates)
    """

    logger = logging.getLogger("OpfLogger")
    logger.info("  generating input e,f values.")

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
