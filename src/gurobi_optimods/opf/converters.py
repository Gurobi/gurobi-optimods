"""
Contains various functions which map data and options from users input to
internal data structures used when formulating models.

- convert_case_to_internal_format: create the 'alldata' struct used internally
  by the formulator from a case data dictionary
- grbmap_coords_from_dict: attach coordinate information (for plotting) to buses
  in the 'alldata' struct
- grbmap_volts_from_dict: attach voltage information (for violation checker) to
  buses in the 'alldata' struct
- build_internal_settings: return internal settings dictionary based on
  user-provided options
"""

import logging
import math

from gurobi_optimods.opf.structs import Branch, Bus, Gen

logger = logging.getLogger(__name__)


def build_internal_settings(
    opftype,
    polar,
    useef,
    usejabr,
    ivtype,
    branchswitching,
    usemipstart,
    minactivebranches,
    useactivelossineqs,
):
    """Take options given to the public API and return a dictionary of settings
    to be used internally when formulating models"""

    # Ensure various defaults are populated.
    settings = {
        "doac": False,
        "dodc": False,
        "doiv": False,
        "dopolar": bool(polar),
        "use_ef": bool(useef),
        "skipjabr": not bool(usejabr),
        "ivtype": "aggressive",
        "branchswitching_mip": False,
        # Formulation for branch-switching where the binary variables simply
        # multiply the continuous variables. Sometimes it works better than
        # branchswitching_mip. Only applicable for AC
        "branchswitching_comp": False,
        "usemipstart": bool(usemipstart),
        "minactivebranches": float(minactivebranches),
        # New linear inequalities developed and implemented by Dan.
        # They are outer approximations of the JABR inequalities
        "useactivelossineqs": bool(useactivelossineqs),
        ### The following settings are internal and not exposed in the API
        # (approximately) fix c, s variables if a voltage solution was read in
        "fixcs": False,
        "fixtolerance": 1.0e-5,
        "usemaxdispersion": False,  # difference between all bus angles is small
        "usemaxphasediff": False,  # difference between 2 adjacent branches is small
        "maxdispersion_deg": 0.0,
        "maxphasediff_deg": 360.0,
        "usequadcostvar": False,
    }

    settings.update(
        {
            "maxdispersion_rad": (math.pi / 180.0) * settings["maxdispersion_deg"],
            "maxphasediff_rad": (math.pi / 180.0) * settings["maxphasediff_deg"],
        }
    )

    # Set the AC/DC/IV model type
    opftype = opftype.lower()
    if opftype in ["ac", "dc", "iv"]:
        settings["do" + opftype] = True
    else:
        raise ValueError(f"Unknown OPF type {opftype}")

    # Sub-type for IV models
    ivtype = ivtype.lower()
    if ivtype in ["plain", "aggressive"]:
        settings["ivtype"] = ivtype
    else:
        raise ValueError(f"Unknown ivtype {ivtype}")

    # Branch switching configurations
    if branchswitching == 0:
        settings["branchswitching_mip"] = False
        settings["branchswitching_comp"] = False
    elif branchswitching == 1:
        settings["branchswitching_mip"] = True
        settings["branchswitching_comp"] = False
    elif branchswitching == 2:
        settings["branchswitching_mip"] = False
        settings["branchswitching_comp"] = True
    else:
        raise ValueError(f"Unknown branchswitching setting {branchswitching}.")

    return settings


def grbmap_coords_from_dict(alldata, coords_dict):
    """
    Maps given case data with given bus coordinates

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param coords_dict: Dictionary holding a mapping of bus index to a given coordinate
    :type coords_dict: dict
    """

    numbuses = alldata["numbuses"]
    buses = alldata["buses"]

    for bnum in range(1, numbuses + 1):
        index = buses[bnum].nodeID
        buses[bnum].lat = coords_dict[index][0]
        buses[bnum].lon = coords_dict[index][1]


def grbmap_volts_from_dict(alldata, volts_dict):
    """
    Maps given case data with given bus input voltages

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param volts_dict: Dictionary holding a mapping of bus index to a given voltage input
    :type volts_dict: dict
    """

    # TODO why not use the same pattern as the coordinate mapper?
    buses = alldata["buses"]
    IDtoCountmap = alldata["IDtoCountmap"]

    for v in volts_dict:
        buscount = IDtoCountmap[v]
        buses[buscount].inputvoltage = 1
        buses[buscount].inputV = volts_dict[v][0]  # Vm
        vadeg = volts_dict[v][1]
        buses[buscount].inputA_rad = vadeg * math.pi / 180.0  # VaRad


def convert_case_to_internal_format(case_dict):
    """
    Converts case format to internal data format, returning an initialised
    'alldata' dictionary. Raises ValueError if there is any invalid data,
    e.g. bus type not allowed, illegal angle for branch, bad bus references.
    """

    # Input data validation
    if len(case_dict["gen"]) != len(case_dict["gencost"]):
        raise ValueError("Invalid input: mismatch between gen and gencost records")

    for gencost in case_dict["gencost"]:
        if gencost["costtype"] != 2:
            raise ValueError("Invalid input: only generator costtype=2 is supported")
        if gencost["n"] != len(gencost["costvector"]):
            raise ValueError(
                "Invalid input: mismatch between gencost.n and costvector length"
            )
        if gencost["n"] > 3:
            raise ValueError(
                "Invalid input: only quadratic and linear cost functions "
                "(gencost.n <= 3) are supported"
            )

    bus_ids = {bus["bus_i"] for bus in case_dict["bus"]}
    if any(branch["fbus"] not in bus_ids for branch in case_dict["branch"]):
        raise ValueError("Unknown bus ID referenced in branch fbus")
    if any(branch["tbus"] not in bus_ids for branch in case_dict["branch"]):
        raise ValueError("Unknown bus ID referenced in branch tbus")
    if any(gen["bus"] not in bus_ids for gen in case_dict["gen"]):
        raise ValueError("Unknown bus ID referenced in generator bus")

    # For each field we create a key value for use internally.
    # Note: index != nodeID (bus_i) for buses.
    # Note: some parts of the code still rely on these indexes being 1..n and
    # the keys being in ascending order in dictionary iterators. Be very careful
    # when trying to change this.
    # Remove isolated buses and their connected branches
    remove_buses = {bus["bus_i"] for bus in case_dict["bus"] if bus["type"] == 4}
    if remove_buses:
        logger.info(f"Removing buses {remove_buses} (bustype=4)")

    # For each field we create a key value for use internally.
    # Note: index != nodeID (bus_i) for buses.
    # Note: some parts of the code still rely on these indexes being 1..n and
    # the keys being in ascending order in dictionary iterators. Be very careful
    # when trying to change this.
    case_dict = {
        "baseMVA": case_dict["baseMVA"],
        "bus": {
            i + 1: dict(bus)
            for i, bus in enumerate(case_dict["bus"])
            if bus["bus_i"] not in remove_buses
        },
        "branch": {
            i + 1: dict(branch)
            for i, branch in enumerate(case_dict["branch"])
            if branch["fbus"] not in remove_buses and branch["tbus"] not in remove_buses
        },
        "gen": {i + 1: dict(gen) for i, gen in enumerate(case_dict["gen"])},
        "gencost": {
            i + 1: dict(gencost) for i, gencost in enumerate(case_dict["gencost"])
        },
        "casename": case_dict["casename"],
    }

    # Manual correction to ratios
    for branch in case_dict["branch"].values():
        if branch["ratio"] == 0.0:
            branch["ratio"] = 1.0

    alldata = {"LP": {}, "MIP": {}}
    alldata["casename"] = case_dict["casename"]
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

        # Could save a lot of pain by using the nodeID as the ... node ID :-)
        # assert nodeID == count
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
                f"Branch # {numbranches} has illegal angle constraints. "
                f"minangle: {minangle} > {maxangle} :maxangle"
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

        # From MATPOWER Generator Data Format
        # https://matpower.org/docs/ref/matpower5.0/caseformat.html
        #
        # 8   status,  >  0 - machine in service
        #             <= 0 - machine out of service
        dgen["status"] == 0 if dgen["status"] <= 0 else 1
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
