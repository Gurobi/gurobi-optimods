import math


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
        # Formulation for branch-switching where the binary variables simply multiply the continuous variables.
        # Sometimes it works better than branchswitching_mip. Only applicable for AC
        "branchswitching_comp": False,
        "usemipstart": bool(usemipstart),
        "minactivebranches": float(minactivebranches),
        # New linear inequalities developed and implemented by Dan.
        # They are outer approximations of the JABR inequalities
        "useactivelossineqs": bool(useactivelossineqs),
        #############################
        # The following settings should currently not be disclosed
        # (TODO: we can't exactly hide these ... what do we mean "not be disclosed")
        # For now keep for us, mainly used for debugging and experimenting with heuristics
        "fixcs": False,  # (approximately) fix c, s variables if a voltage solution was read in
        "fixtolerance": 1.0e-5,
        # Heuristics to help NL solver find a good solution
        "usemaxdispersion": False,  # difference between all bus angles is small
        "usemaxphasediff": False,  # difference between 2 adjacent branches is small
        "maxdispersion_deg": 0.0,
        "maxphasediff_deg": 360.0,
        "usequadcostvar": False,
    }

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
