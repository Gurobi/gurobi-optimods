import math
import cmath
import time
import logging
import scipy
import numpy as np

from .utils import initialize_logger, remove_and_close_handlers


class Bus:
    """
    Describes a bus in a power system, including loads,
    voltage limits, type of bus (generator, reference, etc.), ID and count.
    The data structure also keeps track of branches incident with the bus
    """

    def __init__(
        self,
        count,
        nodeID,
        nodetype,
        Pd,
        Qd,
        Gs,
        Bs,
        area,
        Vm,
        Va,
        Vbase,
        zone,
        Vmax,
        Vmin,
    ):
        self.genidsbycount = []  # array of generator IDs at this bus
        self.frombranchids = {}  # branches where this bus is the 'from' bus
        self.tobranchids = {}  # branches where this bus is the 'to' bus
        self.count = count  # bus count # TODO-Dan Why do we need the count?
        self.nodeID = nodeID  # ID of bus
        self.nodetype = nodetype
        self.Pd = Pd  # active load
        self.Qd = Qd  # reactive load
        self.Gs = Gs  # shunt admittance parameter
        self.Bs = Bs  # shunt admittance parameter
        self.area = area  # we don't use it but keep for consistency with MATPOWER
        self.Vm = Vm  # we don't use it but keep for consistency with MATPOWER
        self.Va = Va  # we don't use it but keep for consistency with MATPOWER
        self.Vbase = Vbase  # voltage base
        self.zone = zone  # we don't use it but keep for consistency with MATPOWER
        self.Vmax = Vmax
        self.Vmin = Vmin
        self.inputvoltage = False
        self.cffvarind = -1  # the following four fields are variable indices
        self.Pinjvarind = -1  # index of Pd solution variable for buses
        self.Qinjvarind = -1  # index of Qd solution variable for buses
        self.vvarind = -1  # index of bus voltage magnitude variable for buses
        self.thetavarind = -1  # index of bus volate angle variable for buses
        self.Pbalance = 0
        self.inputV = 0  # input voltage
        self.inpute = 0  # e-value (input or derived from input voltage solution)
        self.inputf = 0  # e-value (input or derived from input voltage solution)
        self.outdegree = self.indegree = self.degree = 0

        self.lat = -1
        self.lon = -1

    def addgenerator(self, generatorcount, generator):
        logger = logging.getLogger("OpfLogger")
        self.genidsbycount.append(generatorcount)
        loud = False
        if loud:
            logger.info(
                " added generator # "
                + str(generatorcount)
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
    TODO-Dan Please add more description similar to the Bus class
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
    ):
        self.count = count
        self.f = f  # bus ID for from bus
        self.t = t  # bus ID for to bus
        self.count_f = count_f  # count for from bus
        self.count_t = count_t  # count for to bus
        self.r = r  # resistance
        self.x = x  # reactance
        self.bc = bc  # branch charging admittance
        self.count = count  # TODO-Dan Why do we need the count?
        self.rateAmva = rateAmva  # the following three parameters
        self.rateBmva = rateBmva  # describe branch limits
        self.rateCmva = rateCmva
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

    def show(self):  # TODO-Dan This is used nowhere
        logger = logging.getLogger("OpfLogger")
        logger.info(" < " + str(self.f) + " , " + str(self.t) + " > ")
        logger.info(" r " + str(self.r) + " x " + str(self.x) + " bc " + str(self.bc))
        logger.info(" ra " + str(self.ratio) + " ang " + str(self.angle))


class Gen:
    """
    A generator models a complex power injection at a specific bus

    See Matpower manual for more details
    TODO-Dan Could you please at least add a link
    https://matpower.org/docs/MATPOWER-manual.pdf # TODO-Dan is it this one?
    TODO-Dan Please add more description similar to the Bus class if you think it's needed
    """

    def __init__(
        self,
        count,
        nodeID,
        Pg,
        Qg,
        Qmax,
        Qmin,
        Vg,
        mBase,
        status,
        Pmax,
        Pmin,
        Pc1,
        Pc2,
        Qc1min,
        Qc1max,
        Qc2min,
        Qc2max,
        ramp_agc,
        ramp_10,
        ramp_30,
        ramp_q,
        apf,
    ):
        self.count = count  # TODO-Dan Why do we need the count?
        self.nodeID = nodeID  # ID of bus holding gen
        self.Pg = Pg  # active power generation in existing solution
        self.Qg = Qg  # reactive power generation in existing solution
        self.Qmax = Qmax
        self.Qmin = Qmin
        self.Vg = Vg  # not used but kept for completeness
        self.mBase = mBase  # not used but kept for completeness
        self.status = status  # on or off
        self.Pmax = Pmax  # max, min active and reactive power limits
        self.Pmin = Pmin
        self.Pc1 = Pc1  # not used but kept for completeness
        self.Pc2 = Pc2  # not used but kept for completeness
        self.Qc1min = Qc1min  # not used but kept for completeness
        self.Qc1max = Qc1max  # not used but kept for completeness
        self.Qc2min = Qc2min  # not used but kept for completeness
        self.Qc2max = Qc2max  # not used but kept for completeness
        self.ramp_agc = ramp_agc  # not used but kept for completeness
        self.ramp_10 = ramp_10  # not used but kept for completeness
        self.ramp_30 = ramp_30  # not used but kept for completeness
        self.ramp_q = ramp_q  # not used but kept for completeness
        self.apf = apf  # not used but kept for completeness
        self.costlinenum = (
            -1
        )  # TODO-Dan This is used nowhere. Do we need the cost line number?
        self.Pvarind = -1  #
        self.Qvarind = -1

    def addcost(
        self, costvector, linenum
    ):  # TODO-Dan This is not used anywhere. Do we need it?
        self.costvector = costvector
        self.costdegree = len(costvector) - 1
        self.costlinenum = linenum

    def addcost_plain(self, costvector, costtype, startup, shutdown, baseMVA):
        self.costdegree = len(costvector) - 1
        # Re-scale costs
        for j in range(self.costdegree + 1):
            costvector[j] *= baseMVA ** (self.costdegree - j)
        self.costvector = costvector
        self.costtype = costtype
        self.startup = startup
        self.shutdown = shutdown
        self.costlinenum = -1

    def showcostvector(self):  # TODO-Dan This is used nowhere
        logger = logging.getLogger("OpfLogger")
        logger.info(self.costvector)
        for i in range(0, self.costdegree + 1):
            logger.info(i, self.costvector[i], " ", end="")


def read_case(alldata, case_dict):
    """
    Fills alldata dictionary out of a user-given case dictionary

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    case_dict : dictionary
        Dictionary holding all case relevant data
    """

    logger = logging.getLogger("OpfLogger")
    buses = {}
    IDtoCountmap = {}
    slackbus = -1
    numbuses = 0
    numPload = 0
    sumPd = sumQd = 0
    numisolated = 0

    logger.info("Building case data structures from dictionary.")

    dict_buses = case_dict["bus"]
    baseMVA = alldata["baseMVA"] = case_dict["baseMVA"]

    logger.info("Buses.")

    for dbus in dict_buses.values():
        numbuses += 1
        count = numbuses
        nodeID = dbus["bus_i"]
        nodetype = dbus["type"]

        if nodetype != 1 and nodetype != 2 and nodetype != 3 and nodetype != 4:
            raise ValueError("Bad bus %s has type %s." % (count, nodetype))

        if nodetype == 3:
            slackbus = count
            logger.info("    Bus %d ID %d is the reference bus." % (count, nodeID))
            alldata["refbus"] = count

        if nodetype == 4:
            numisolated += 1

        Pd = dbus["Pd"]
        Qd = dbus["Qd"]
        Gs = dbus["Gs"]
        Bs = dbus["Bs"]
        area = dbus["area"]
        Vm = dbus["Vm"]
        Va = dbus["Va"]
        Vbase = dbus["baseKV"]
        zone = dbus["zone"]
        Vmax = dbus["Vmax"]
        Vmin = dbus["Vmin"]

        buses[numbuses] = Bus(
            numbuses,
            nodeID,
            nodetype,
            Pd / baseMVA,
            Qd / baseMVA,
            Gs / baseMVA,
            Bs / baseMVA,
            area,
            Vm,
            Va,
            Vbase,
            zone,
            Vmax,
            Vmin,
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

    dict_branches = case_dict["branch"]
    branches = {}

    for dbranch in dict_branches.values():
        # print(dbranch)
        numbranches += 1
        brcnt1 = numbranches
        f = dbranch["fbus"]
        t = dbranch["tbus"]
        r = dbranch["r"]
        x = dbranch["x"]
        bc = dbranch["b"]
        rateA = dbranch["rateA"]
        rateB = dbranch["rateB"]
        rateC = dbranch["rateC"]
        ratio = dbranch["ratio"]
        if ratio == 0:
            ratio = 1.0  # to be sure
        angle = dbranch["angle"]
        status = dbranch["status"]
        minangle = dbranch["angmin"]
        maxangle = dbranch["angmax"]
        if maxangle < minangle:
            raise ValueError("Branch # %d has illegal angle constraints." % numbranches)

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
    dict_gens = case_dict["gen"]
    summaxgenP = summaxgenQ = 0
    numgens = 0

    for dgen in dict_gens.values():
        numgens += 1
        gencount1 = numgens
        nodeID = dgen["bus"]
        Pg = dgen["Pg"]
        Qg = dgen["Qg"]
        Qmax = dgen["Qmax"]
        Qmin = dgen["Qmin"]
        Vg = dgen["Vg"]
        mBase = dgen["mBase"]
        status = dgen["status"]
        Pmax = dgen["Pmax"]
        Pmin = dgen["Pmin"]
        Pc1 = dgen["Pc1"]
        Pc2 = dgen["Pc2"]
        Qc1min = dgen["Qc1min"]
        Qc1max = dgen["Qc1max"]
        Qc2min = dgen["Qc2min"]
        Qc2max = dgen["Qc2max"]
        ramp_agc = dgen["ramp_agc"]
        ramp_10 = dgen["ramp_10"]
        ramp_30 = dgen["ramp_30"]
        ramp_q = dgen["ramp_q"]
        apf = dgen["apf"]

        idgencount1 = -1
        if nodeID in IDtoCountmap.keys():
            idgencount1 = IDtoCountmap[nodeID]
            gens[gencount1] = Gen(
                gencount1,
                nodeID,
                Pg,
                Qg,
                Qmax / baseMVA,
                Qmin / baseMVA,
                Vg,
                mBase,
                status,
                Pmax / baseMVA,
                Pmin / baseMVA,
                Pc1,
                Pc2,
                Qc1min,
                Qc1max,
                Qc2min,
                Qc2max,
                ramp_agc,
                ramp_10,
                ramp_30,
                ramp_q,
                apf,
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
    gencoststruct = case_dict["gencost"]

    for count in range(1, alldata["numgens"] + 1):
        gens[count].addcost_plain(
            gencoststruct[count]["costvector"],
            gencoststruct[count]["costtype"],
            gencoststruct[count]["startup"],
            gencoststruct[count]["shutdown"],
            baseMVA,
        )


def read_case_file(casefile):
    """
    Reads case data from an .m file following MATPOWER standards and
    returns an OptiMod compatible dictionary holding all relevant case data

    Parameters
    ----------
    casefile : string
        Path to .m file holding all case relevant information

    Returns
    -------
    dictionary
        Dictionary holding all case relevant data
    """

    logger, handlers = initialize_logger("CaseReadingLogger")

    logger.info("Reading case file %s." % casefile)
    starttime = time.time()
    f = open(casefile, "r")
    lines = f.readlines()
    f.close()

    # Read all lines of the casefile
    logger.info("Building case dictionary.")
    case_dict = read_case_build_dict(lines)
    endtime = time.time()
    logger.info("Reading and building time: %f s." % (endtime - starttime))
    logger.info("")

    remove_and_close_handlers(logger, handlers)

    return case_dict


def read_case_build_dict(lines):
    """
    Reads thru all lines of a previously read-in case file in .m format and
    returns an OptiMod compatible dictionary holding all relevant case data
    This is a helper function for translating casefiles into dict format

    Parameters
    ----------
    lines : list
        List holding all lines of a previously read-in case file

    Returns
    -------
    dictionary
        Dictionary holding all case relevant data
    """

    logger = logging.getLogger("CaseReadingLogger")
    numlines = len(lines)
    lookingforbus = 1
    linenum = 2
    baseMVA = 0.0

    case_dict = {}

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
            case_dict["bus"] = buses
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

                    Pd = float(thisline[2])
                    Qd = float(thisline[3])
                    Gs = float(thisline[4])
                    Bs = float(thisline[5])
                    area = float(thisline[6])
                    Vm = float(thisline[7])
                    Va = float(thisline[8])
                    Vbase = float(thisline[9])
                    zone = float(thisline[10])
                    Vmax = float(thisline[11])
                    # Trim ; if present
                    Vmin = thisline[12]
                    if Vmin[len(Vmin) - 1] == ";":
                        Vmin = Vmin[: len(Vmin) - 1]
                    Vmin = float(Vmin)

                    if nodetype == 3:
                        logger.info(
                            "    Bus %d ID %d is the reference bus."
                            % (numbuses, nodeID)
                        )

                else:
                    nodeID = int(thisline[1])
                    nodetype = int(thisline[2])
                    logger.info(
                        "bus %d nodeID %d is isolated and has type %d."
                        % (numbuses, nodeID, nodetype)
                    )
                    logger.info("   setting it to type 4.")
                    nodetype = 4

                    Pd = Qd = Gs = Bs = area = Vm = Va = 0
                    Vbase = float(thisline[10])
                    zone = Vmax = Vmin = 0

                buses[numbuses] = {}
                buses[numbuses]["bus_i"] = nodeID
                buses[numbuses]["type"] = nodetype
                buses[numbuses]["Pd"] = Pd
                buses[numbuses]["Qd"] = Qd
                buses[numbuses]["Gs"] = Gs
                buses[numbuses]["Bs"] = Bs
                buses[numbuses]["area"] = area
                buses[numbuses]["Vm"] = Vm
                buses[numbuses]["Va"] = Va
                buses[numbuses]["baseKV"] = Vbase
                buses[numbuses]["zone"] = zone
                buses[numbuses]["Vmax"] = Vmax
                buses[numbuses]["Vmin"] = Vmin

                linenum += 1

            if lookingforendofbus:
                raise ValueError("Could not find bus data section.")

        elif theword == "mpc.gen":
            case_dict["gen"] = gens = {}
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
                Qmax = float(thisline[3])
                Qmin = float(thisline[4])
                Vg = float(thisline[5])
                mBase = float(thisline[6])
                status = int(thisline[7])
                Pmax = float(thisline[8])
                Pmin = float(thisline[9])
                Pc1 = float(thisline[10])
                Pc2 = float(thisline[11])
                Qc1min = float(thisline[12])
                Qc1max = float(thisline[13])
                Qc2min = float(thisline[14])
                Qc2max = float(thisline[15])
                ramp_agc = float(thisline[16])
                ramp_10 = float(thisline[17])
                ramp_30 = float(thisline[18])
                ramp_q = float(thisline[19])
                # Trim ; if present
                apf = thisline[20]
                if apf[len(apf) - 1] == ";":
                    apf = apf[: len(apf) - 1]
                apf = float(apf)

                if status <= 0:
                    status = 0
                else:
                    status = 1

                gens[gencount1] = {}
                gens[gencount1]["bus"] = nodeID
                gens[gencount1]["Pg"] = Pg
                gens[gencount1]["Qg"] = Qg
                gens[gencount1]["Qmax"] = Qmax
                gens[gencount1]["Qmin"] = Qmin
                gens[gencount1]["Vg"] = Vg
                gens[gencount1]["mBase"] = mBase
                gens[gencount1]["status"] = status
                gens[gencount1]["Pmax"] = Pmax
                gens[gencount1]["Pmin"] = Pmin
                gens[gencount1]["Pc1"] = Pc1
                gens[gencount1]["Pc2"] = Pc2
                gens[gencount1]["Qc1min"] = Qc1min
                gens[gencount1]["Qc1max"] = Qc1max
                gens[gencount1]["Qc2min"] = Qc2min
                gens[gencount1]["Qc2max"] = Qc2max
                gens[gencount1]["ramp_agc"] = ramp_agc
                gens[gencount1]["ramp_10"] = ramp_10
                gens[gencount1]["ramp_30"] = ramp_30
                gens[gencount1]["ramp_q"] = ramp_q
                gens[gencount1]["apf"] = apf

                linenum += 1
            # Finished reading gen section
            if lookingforendofgen:
                raise ValueError("Could not find end of generator section.")

        elif theword == "mpc.branch":
            case_dict["branch"] = branches = {}
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
                branches[numbranches]["fbus"] = f
                branches[numbranches]["tbus"] = t
                branches[numbranches]["r"] = r
                branches[numbranches]["x"] = x
                branches[numbranches]["b"] = bc
                branches[numbranches]["rateA"] = rateA
                branches[numbranches]["rateB"] = rateB
                branches[numbranches]["rateC"] = rateC
                branches[numbranches]["ratio"] = ratio
                branches[numbranches]["angle"] = angle
                branches[numbranches]["status"] = status
                branches[numbranches]["angmin"] = minangle
                branches[numbranches]["angmax"] = maxangle

                linenum += 1

            # Finished reading branch section
            if lookingforendofbranch:
                raise ValueError("Could not find end of branch section.")

        elif theword == "mpc.gencost":
            case_dict["gencost"] = gencoststruct = {}

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
                    startup = thisline[1]
                    shutdown = thisline[2]
                    degree = int(thisline[3])

                    gencoststruct[gencostcount]["costtype"] = costtype
                    gencoststruct[gencostcount]["n"] = degree
                    gencoststruct[gencostcount]["startup"] = startup
                    gencoststruct[gencostcount]["shutdown"] = shutdown

                    costvector = [float(thisline[j]) for j in range(4, 4 + degree - 1)]
                    # Trim ; at the end of line
                    lastcost = thisline[-1]
                    if lastcost[len(lastcost) - 1] == ";":
                        lastcost = lastcost[: len(lastcost) - 1]
                    lastcost = float(lastcost)
                    costvector.append(lastcost)
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


def read_case_file_mat(casefile):
    """
    Reads case data from an .mat file following MATPOWER standards and
    returns an OptiMod compatible dictionary holding all relevant case data

    Parameters
    ----------
    casefile : string
        Path to .mat file holding all case relevant information

    Returns
    -------
    dictionary
        A dictionary holding all case relevant data
    """

    logger, handlers = initialize_logger("CaseReadingLogger")
    starttime = time.time()
    case_dict = {}
    case_dict["refbus"] = -1

    logger.info("Reading case file %s." % casefile)
    mat = loadmat(casefile)

    if "mpc" not in mat.keys():
        raise ValueError("Provided .mat file does not have an mpc field")

    mpc = mat["mpc"]

    for x in ["baseMVA", "bus", "gen", "branch", "gencost"]:
        if x not in mpc.keys():
            raise ValueError("Provided .mat file does not have a %s field" % x)

    case_dict["baseMVA"] = mpc["baseMVA"]
    mpcbuses = mpc["bus"]
    mpcgen = mpc["gen"]
    mpcbranch = mpc["branch"]
    mpcgencost = mpc["gencost"]

    print(mat)
    numbuses = 0
    buses = {}
    for b in mpcbuses:
        numbuses += 1
        buses[numbuses] = {}
        buses[numbuses]["bus_i"] = int(b[0])
        buses[numbuses]["type"] = int(b[1])
        buses[numbuses]["Pd"] = b[2]
        buses[numbuses]["Qd"] = b[3]
        buses[numbuses]["Gs"] = b[4]
        buses[numbuses]["Bs"] = b[5]
        buses[numbuses]["area"] = b[6]
        buses[numbuses]["Vm"] = b[7]
        buses[numbuses]["Va"] = b[8]
        buses[numbuses]["baseKV"] = b[9]
        buses[numbuses]["zone"] = b[10]
        buses[numbuses]["Vmax"] = b[11]
        buses[numbuses]["Vmin"] = b[12]
        if buses[numbuses]["type"] == 3:
            slackbus = buses[numbuses]["bus_i"]
            refbus = numbuses
            logger.info("    Slack bus: %d" % slackbus)
            logger.info(
                "    Bus %d ID %d is the reference bus."
                % (numbuses, buses[numbuses]["bus_i"])
            )

    case_dict["bus"] = buses

    numgens = 0
    gens = {}
    for g in mpcgen:
        numgens += 1
        gens[numgens] = {}
        gens[numgens]["bus"] = int(g[0])
        gens[numgens]["Pg"] = g[1]
        gens[numgens]["Qg"] = g[2]
        gens[numgens]["Qmax"] = g[3]
        gens[numgens]["Qmin"] = g[4]
        gens[numgens]["Vg"] = g[5]
        gens[numgens]["mBase"] = g[6]
        gens[numgens]["status"] = 0 if g[7] <= 0 else 1
        gens[numgens]["Pmax"] = g[8]
        gens[numgens]["Pmin"] = g[9]
        gens[numgens]["Pc1"] = g[10]
        gens[numgens]["Pc2"] = g[11]
        gens[numgens]["Qc1min"] = g[12]
        gens[numgens]["Qc1max"] = g[13]
        gens[numgens]["Qc2min"] = g[14]
        gens[numgens]["Qc2max"] = g[15]
        gens[numgens]["ramp_agc"] = g[16]
        gens[numgens]["ramp_10"] = g[17]
        gens[numgens]["ramp_30"] = g[18]
        gens[numgens]["ramp_q"] = g[19]
        gens[numgens]["apf"] = g[20]

    case_dict["gen"] = gens

    numbranches = 0
    branches = {}
    for b in mpcbranch:
        numbranches += 1
        branches[numbranches] = {}
        branches[numbranches]["fbus"] = b[0]
        branches[numbranches]["tbus"] = b[1]
        branches[numbranches]["r"] = b[2]
        branches[numbranches]["x"] = b[3]
        branches[numbranches]["b"] = b[4]
        branches[numbranches]["rateA"] = b[5]
        branches[numbranches]["rateB"] = b[6]
        branches[numbranches]["rateC"] = b[7]
        branches[numbranches]["ratio"] = 1.0 if b[8] == 0.0 else b[8]
        branches[numbranches]["angle"] = b[9]
        branches[numbranches]["status"] = b[10]
        branches[numbranches]["angmin"] = b[11]
        branches[numbranches]["angmax"] = b[12]

    case_dict["branch"] = branches

    numgencosts = 0
    gencosts = {}
    for g in mpcgencost:
        numgencosts += 1
        gencosts[numgencosts] = {}
        gencosts[numgencosts]["costtype"] = g[0]
        gencosts[numgencosts]["startup"] = g[1]
        gencosts[numgencosts]["shutdown"] = g[2]
        degree = int(g[3])
        gencosts[numgencosts]["n"] = degree
        costvector = [g[j] for j in range(4, 4 + degree)]
        gencosts[numgencosts]["costvector"] = costvector

    if numgencosts > numgens:
        raise ValueError(
            "Read %d gen costs but only %d generators." % (numgencosts, numgens)
        )

    case_dict["gencost"] = gencosts

    endtime = time.time()
    logger.info("Reading and building time: %f s." % (endtime - starttime))
    logger.info("")
    remove_and_close_handlers(logger, handlers)

    return case_dict


def loadmat(filename):
    """
    This function should be called instead of direct scipy.io.loadmat
    as it cures the problem of not properly recovering Python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects

    Parameters
    ----------
    filename : string
        Path to .mat file holding all case relevant information

    Returns
    -------
    dictionary
        Dictionary holding all case relevant data read from the given .mat file
    """
    data = scipy.io.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)


def _check_keys(dict):
    """
    Checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries

    Parameters
    ----------
    dict : dictionary
        Dictionary holding all case data read from a .mat file by the
        scipy.io.loadmat function

    Returns
    -------
    dictionary
        Dictionary holding all case relevant data read from the given .mat file
    """
    for key in dict:
        if isinstance(dict[key], scipy.io.matlab.mat_struct):
            dict[key] = _todict(dict[key])
    return dict


def _todict(matobj):
    """
    A recursive function which constructs from matobjects nested dictionaries

    Parameters
    ----------
    matobj : scipy.io.matlab.mat_struct
        Placeholder for holding read data from structs

    Returns
    -------
    dictionary
        A possibly nested dictionary with all case relevant data
    """
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, scipy.io.matlab.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict


def turn_opf_dict_into_mat_file(solution, filename):
    """
    Writes a .mat file out of and OPF solution dictionary

    Parameters
    ----------
    solution : dictionary
        OPF solution dictionary
    filename : string
        Name of .mat file where to write the solution data
    """

    # Buses
    buses = solution["bus"]
    matrix = []
    for bus in buses.values():
        matrix.append(list(bus.values()))

    solution["bus"] = np.array(matrix)
    # Generators
    gens = solution["gen"]
    matrix = []
    for gen in gens.values():
        matrix.append(list(gen.values()))

    solution["gen"] = np.array(matrix)
    # Branches
    branches = solution["branch"]
    matrix = []
    for branch in branches.values():
        matrix.append(list(branch.values()))

    solution["branch"] = np.array(matrix)
    # Generator costs
    gencosts = solution["gencost"]
    matrix = []
    # For generator costs, the last dictionary entry is a list of values
    for genc in gencosts.values():
        l = list(genc.values())
        costvector = l[-1]
        l.pop()
        for item in costvector:
            l.append(item)

        matrix.append(l)

    solution["gencost"] = np.array(matrix)
    # Write mat file
    scipy.io.savemat(filename, {"result": solution})
