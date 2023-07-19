import cmath
import logging
import math

logger = logging.getLogger(__name__)


class Bus:
    """
    Describes a bus in a power system, including loads,
    voltage limits, type of bus (generator, reference, etc.), ID and count.
    The data structure also keeps track of branches incident with the bus

    :param count: Number of bus in the order of reading it in
    :type count: int
    :param nodeID: Bus ID
    :type nodeID: int
    :param nodetype: Bus type: PQ bus = 1, PV bus = 2, reference bus = 3, isolated bus = 4
    :type nodetype: int
    :param Pd: Real power demand
    :type Pd: float
    :param Qd: reactive power demand
    :type Qd: float
    :param Gs: shunt conductance
    :type Gs: float
    :param Bs: shunt susceptance
    :type Bs: float
    :param area: area number (currently not used)
    :type area: int
    :param Vm: Voltage magnitude
    :type Vm: float
    :param Va: Voltage angle
    :type Va: float
    :param Vbase: Base voltage
    :type Vbase: float
    :param zone: loss zone (currently not used)
    :type zone: int
    :param Vmax: Maximum voltage magnitude
    :type Vmax: float
    :param Vmin: Minimum voltage magnitude
    :type Vmin: float
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
        """Constructor method"""
        self.genidsbycount = []  # array of generator IDs at this bus
        self.frombranchids = {}  # branches where this bus is the 'from' bus
        self.tobranchids = {}  # branches where this bus is the 'to' bus
        self.count = count  # bus count, needed to identify bus without ID
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
        self.Vmax = Vmax  # max voltage
        self.Vmin = Vmin  # min voltage
        self.inputvoltage = False
        self.inputV = 0  # input voltage
        self.inpute = 0  # e-value (input or derived from input voltage solution)
        self.inputf = 0  # e-value (input or derived from input voltage solution)
        self.outdegree = self.indegree = self.degree = 0
        # Coordinates for plotting
        self.lat = -1
        self.lon = -1

    def addgenerator(self, generatorcount, generator):
        """
        Adds a generator to the bus

        :param generatorcount: Number of generator added to the bus
        :type generatorcount: int
        :parameter generator: Generator to be added to the bus
        :type generator: Gen
        """
        self.genidsbycount.append(generatorcount)

    def addfrombranch(self, id):
        """
        Adds a "from" branch to the bus

        :param id: ID of "from" branch
        :type id: int
        """
        quant = len(self.frombranchids)
        self.frombranchids[quant] = id
        self.outdegree += 1
        self.degree += 1

    def addtobranch(self, id):
        """
        Adds a "to" branch to the bus

        :param id: ID of "to" branch
        :type id: int
        """
        quant = len(self.tobranchids)
        self.tobranchids[quant] = id
        self.indegree += 1
        self.degree += 1


class Branch:
    """
    Describes a branch in a power flow network which connects two buses.

    :param count: Number of branch in the order of reading it in
    :type count: int
    :param f: "from" bus ID
    :type f: int
    :param count_f: count attribute value of "from" bus
    :type count_f: int
    :param t: "to" bus ID
    :type t: int
    :param count_t: count attribute value of "to" bus
    :type count_t: int
    :param r: resistance
    :type r: float
    :param x: reactance
    :type x: float
    :param bc: total line charging susceptance
    :type bc: float
    :param rateAmva: MVA long term rating
    :type rateAmva: float
    :param rateBmva: MVA short term rating
    :type rateBmva: float
    :param rateCmva: MVA emergency term rating
    :type rateCmva: float
    :param ratio: tranformer off nominal turns ratio
    :type ratio: int
    :param angle: transformer phase shift angle (degrees)
    :type angle: float
    :param maxangle: maximum angle difference (degrees)
    :param status: States whether branch is on or off
    :type status: int
    :type maxangle: float
    :param minangle: minimum angle difference (degrees)
    :type minangle: float
    :param defaultlimit: default limit on branch voltage
    :type defaultlimit: float
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
        status,
        maxangle,
        minangle,
        defaultlimit,
    ):
        """Constructor method"""
        self.count = count
        self.f = f  # bus ID for from bus
        self.t = t  # bus ID for to bus
        self.count_f = count_f  # count for from bus
        self.count_t = count_t  # count for to bus
        self.r = r  # resistance
        self.x = x  # reactance
        self.bc = bc  # branch charging admittance
        self.count = count  # needed to identify branch
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

        if maxangle == 360 or maxangle == 0:
            self.maxangle_rad = 2 * math.pi

        if minangle == -360 or minangle == 0:
            self.minangle_rad = -2 * math.pi

        self.invratio2 = invratio2 = 1 / ratio**2
        self.multtf = multtf = 1 / (ratio * cmath.exp(1j * self.angle_rad))
        self.multft = multft = 1 / (ratio * cmath.exp(-1j * self.angle_rad))
        # print 'multtf', multtf
        self.status = status
        self.z = z = r + x * 1j
        self.y = y = 1 / z
        self.ynorm2 = (y.real) * (y.real) + (y.imag) * (y.imag)
        self.Yff = (y + bc / 2 * 1j) * invratio2
        self.Yft = -y * multft
        self.Ytf = -y * multtf
        self.Ytt = y + bc / 2 * 1j
        self.Gff = (self.Yff).real
        self.Bff = (self.Yff).imag
        self.Gft = (self.Yft).real
        self.Bft = (self.Yft).imag
        self.Gtf = (self.Ytf).real
        self.Btf = (self.Ytf).imag
        self.Gtt = (self.Ytt).real
        self.Btt = (self.Ytt).imag

        self.isacline = (ratio == 1) and (
            self.angle_rad == 0
        )  # Currently not used, keep for future work
        self.nongaining = (
            (self.Gff >= 0)
            and (self.Gtt >= 0)
            and (self.Gff >= -self.Gft)
            and (self.Gtt >= -self.Gtf)
        )

        self.inputcs = False  # Parameter for future work
        self.inputc = 2  # Parameter for future work
        self.inputs = 2  # Parameter for future work

    def show(self):
        """
        Prints a branch. Could be useful for debugging
        """
        logger.info(f" < {self.f} , {self.t} > ")
        logger.info(f" r {self.r} x {self.x} bc {self.bc}")
        logger.info(f" ra {self.ratio} ang {self.angle}")


class Gen:
    """
    A generator models a complex power injection at a specific bus

    See Matpower manual for more details
    https://matpower.org/docs/MATPOWER-manual.pdf

    :param count: Number of generator in the order of reading it in
    :type count: int
    :param nodeID: Generator ID
    :type nodeID: int
    :param Pg: Real power output
    :type Pg: float
    :param Qg: Reactive power output
    :type Qg: float
    :param Qmax: Maximum reactive power output
    :type Qmax: float
    :param Qmin: Minimum reactive power output
    :type Qmin: float
    :param Vg: Voltage magnitude setpoint
    :type Vg: float
    :param mBase: Total MVA base of this machine
    :type mBase: float
    :param status: > 0 generator in service, <= 0 out of service
    :type status: int
    :param Pmax: Maximum real power output
    :type Pmax: float
    :param Pmin: Minimum real power output
    :type Pmin: float
    :param Pc1: Lower real power output of PQ capability curve
    :type Pc1: float
    :param Pc2: Upper real power output of PQ capability curve
    :type Pc2: float
    :param Qc1min: Minimum reactive power output at Pc1
    :type Qc1min: float
    :param Qc1max: Maximum reactive power output at Pc1
    :type Qc1max: float
    :param Qc2min: Minimum reactive power output at Pc2
    :type Qc2min: float
    :param Qc2max: Maximum reactive power output at Pc2
    :type Qc2max: float
    :param ramp_agc: Ramp rate for load following/AGC
    :type ramp_agc: float
    :param ramp_10: Ramp rate for 10 minute reserves
    :type ramp_10: float
    :param ramp_30: Ramp rate for 30 minute reserves
    :type ramp_30: float
    :param ramp_q: Ramp rate for reactive power
    :type ramp_q: float
    :param apf: Area participation factor
    :type apf: float
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
        """Constructor method"""
        self.count = count  # needed to identify gen
        self.nodeID = nodeID  # ID of bus holding gen
        self.Pg = Pg  # active power generation in existing solution
        self.Qg = Qg  # reactive power generation in existing solution
        self.Qmax = Qmax  # max, min reactive power limits
        self.Qmin = Qmin
        self.Vg = Vg  # not used but kept for completeness
        self.mBase = mBase  # not used but kept for completeness
        self.status = status  # on or off
        self.Pmax = Pmax  # max, min active power limits
        self.Pmin = Pmin
        self.Pc1 = Pc1  # not used but kept for completeness with MATPOWER
        self.Pc2 = Pc2  # not used but kept for completeness with MATPOWER
        self.Qc1min = Qc1min  # not used but kept for completeness with MATPOWER
        self.Qc1max = Qc1max  # not used but kept for completeness with MATPOWER
        self.Qc2min = Qc2min  # not used but kept for completeness with MATPOWER
        self.Qc2max = Qc2max  # not used but kept for completeness with MATPOWER
        self.ramp_agc = ramp_agc  # not used but kept for completeness with MATPOWER
        self.ramp_10 = ramp_10  # not used but kept for completeness with MATPOWER
        self.ramp_30 = ramp_30  # not used but kept for completeness with MATPOWER
        self.ramp_q = ramp_q  # not used but kept for completeness with MATPOWER
        self.apf = apf  # not used but kept for completeness with MATPOWER

    def addcost_plain(self, costvector, costtype, startup, shutdown, baseMVA):
        """
        Updates generator cost data.
        Scales the values by given degree and baseMVA

        :param costvector: List of generator costs
        :type costvector: list
        :param costtype: Generator cost type
        :type costtype: int
        :param startup: Startup cost in US dollar
        :type startup: float
        :param shutdown: Shutdown cost in US dollar
        :type shutdown: float
        :param baseMVA: base voltage for scaling
        :type baseMVA: float
        """
        self.costdegree = len(costvector) - 1
        # Re-scale costs
        for j in range(self.costdegree + 1):
            costvector[j] *= baseMVA ** (self.costdegree - j)
        self.costvector = costvector
        self.costtype = costtype
        self.startup = startup
        self.shutdown = shutdown

    def showcostvector(self):
        """
        Prints the generator cost vector. Could be useful for debugging.
        """
        logger.info(self.costvector)
        for i in range(0, self.costdegree + 1):
            logger.info(i, self.costvector[i], " ", end="")


def convert_case_to_internal_format(case_dict):
    """
    Converts case format to internal data format, returning an initialised
    'alldata' dictionary. Raises ValueError if there is any invalid data,
    e.g. bus type not allowed, illegal angle for branch, bad bus references.
    """

    # Pre-conversion needed for hacky approach to indexing
    import copy

    case_dict = copy.deepcopy(case_dict)

    def matlabify(values):
        # 0..n-1 list -> 1..n dict
        return dict(zip(range(1, len(values) + 1), values))

    for branch in case_dict["branch"]:
        if branch["ratio"] == 0.0:
            branch["ratio"] = 1.0
    for field in ["bus", "gen", "branch", "gencost"]:
        case_dict[field] = matlabify(case_dict[field])
    # case_dict is now matlab-like

    alldata = {"LP": {}, "MIP": {}}

    logger.info("Building case data structures from dictionary.")
    baseMVA = alldata["baseMVA"] = case_dict["baseMVA"]
    # Buses
    logger.info("Buses.")
    buses = {}
    IDtoCountmap = {}
    slackbus = refbus = -1
    numbuses = 0
    numPload = 0
    sumPd = sumQd = 0
    numisolated = 0
    dict_buses = case_dict["bus"]

    for dbus in dict_buses.values():
        numbuses += 1
        count = numbuses
        nodeID = dbus["bus_i"]
        nodetype = dbus["type"]
        Pd = dbus["Pd"]
        Qd = dbus["Qd"]

        if nodetype not in [1, 2, 3, 4]:
            raise ValueError(
                f"Bus {count} has type {nodetype}. Only bus types [1,2,3,4] allowed."
            )

        if nodetype == 3:
            slackbus = refbus = count
            logger.info(f"    Bus {count} ID {nodeID} is the reference bus.")

        if nodetype == 4:
            numisolated += 1

        if nodetype in [1, 2, 3]:
            sumPd += Pd
            sumQd += Qd

        buses[numbuses] = Bus(
            numbuses,
            nodeID,
            nodetype,
            Pd / baseMVA,
            Qd / baseMVA,
            dbus["Gs"] / baseMVA,
            dbus["Bs"] / baseMVA,
            dbus["area"],
            dbus["Vm"],
            dbus["Va"],
            dbus["baseKV"],
            dbus["zone"],
            dbus["Vmax"],
            dbus["Vmin"],
        )

        IDtoCountmap[nodeID] = count
        numPload += Pd > 0

    alldata["buses"] = buses
    alldata["numbuses"] = numbuses
    alldata["sumPd"] = sumPd
    alldata["sumQd"] = sumQd
    alldata["IDtoCountmap"] = IDtoCountmap
    alldata["slackbus"] = slackbus
    alldata["refbus"] = refbus

    logger.info(f"    sumloadPd {sumPd} numPload {numPload}")
    logger.info(f"    sumloadQd {sumQd}")

    if slackbus < 0:
        logger.info("    Could not find slack bus.")

    logger.info(f"    {numbuses} buses")
    if numisolated > 0:
        logger.info(f"    Isolated buses: {numisolated}")

    # Branches
    logger.info("Branches.")
    branches = {}
    defaultlimit = 1e20
    numbranches = 0
    activebranches = 0
    dict_branches = case_dict["branch"]

    for dbranch in dict_branches.values():
        # print(dbranch)
        numbranches += 1
        brcnt1 = numbranches
        ratio = dbranch["ratio"]
        minangle = dbranch["angmin"]
        maxangle = dbranch["angmax"]
        f = dbranch["fbus"]
        t = dbranch["tbus"]
        count_f = IDtoCountmap[f]
        count_t = IDtoCountmap[t]

        # TODO why do we adjust this manually? Is this invalid input data?
        if ratio == 0:
            ratio = 1.0  # to be sure

        if maxangle < minangle:
            raise ValueError(
                f"Branch # {numbranches} has illegal angle constraints. minangle: {minangle} > {maxangle} :maxangle"
            )

        branches[numbranches] = Branch(
            brcnt1,
            f,
            count_f,
            t,
            count_t,
            dbranch["r"],
            dbranch["x"],
            dbranch["b"],
            dbranch["rateA"] / baseMVA,
            dbranch["rateB"] / baseMVA,
            dbranch["rateC"] / baseMVA,
            ratio,
            dbranch["angle"],
            dbranch["status"],
            maxangle,
            minangle,
            defaultlimit,
        )
        activebranches += 1
        buses[count_f].addfrombranch(brcnt1)
        buses[count_t].addtobranch(brcnt1)

    alldata["branches"] = branches
    alldata["numbranches"] = numbranches
    logger.info(f"    Numbranches: {numbranches} active: {activebranches}")

    # Generators
    logger.info("Generators.")
    gens = {}
    summaxgenP = summaxgenQ = 0
    numgens = 0
    dict_gens = case_dict["gen"]

    for dgen in dict_gens.values():
        numgens += 1
        gencount1 = numgens
        nodeID = dgen["bus"]
        Qmax = dgen["Qmax"]
        Pmax = dgen["Pmax"]
        idgencount1 = IDtoCountmap[nodeID]

        if nodeID not in IDtoCountmap.keys():
            raise ValueError(f"Generator # {gencount1} in nonexistent bus ID {nodeID}.")

        if buses[idgencount1].nodetype in [2, 3]:  # But not 4
            summaxgenP += Pmax
            summaxgenQ += Qmax

        # TODO status > 0 generator in service, <= 0 out of service
        # Is this how the spec works? The old reader code did this:
        #
        # gens[numgens]["status"] = 0 if g[7] <= 0 else 1
        #
        # but we don't have any tests that hit it.
        assert dgen["status"] in (0, 1)

        gens[gencount1] = Gen(
            gencount1,
            nodeID,
            dgen["Pg"],
            dgen["Qg"],
            Qmax / baseMVA,
            dgen["Qmin"] / baseMVA,
            dgen["Vg"],
            dgen["mBase"],
            dgen["status"],
            Pmax / baseMVA,
            dgen["Pmin"] / baseMVA,
            dgen["Pc1"],
            dgen["Pc2"],
            dgen["Qc1min"],
            dgen["Qc1max"],
            dgen["Qc2min"],
            dgen["Qc2max"],
            dgen["ramp_agc"],
            dgen["ramp_10"],
            dgen["ramp_30"],
            dgen["ramp_q"],
            dgen["apf"],
        )

        buses[idgencount1].addgenerator(gencount1, gens[gencount1])

    alldata["gens"] = gens
    alldata["numgens"] = len(gens)
    busgencount = 0

    for bus in buses.values():
        busgencount += len(bus.genidsbycount) > 0

    alldata["busgencount"] = busgencount
    alldata["summaxgenP"] = summaxgenP
    alldata["summaxgenQ"] = summaxgenQ

    logger.info(f"    Number of generators: {alldata['numgens']}")
    logger.info(f"    Number of buses with gens: {busgencount}")
    logger.info(f"    summaxPg {summaxgenP} summaxQg {summaxgenQ}")

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

    return alldata
