from enum import Enum


class OpfType(Enum):
    AC = "AC"
    DC = "DC"
    IV = "IV"


def get_default_optimization_settings(casefile):
    """Returns a dictionary holding all default settings"""

    settings = {
        "casefilename": None,
        "voltsfilename": None,  # TODO-Dan could you provide an example of how to use this?
        "lpfilename": "grbopf.lp",  # TODO default should be None
        "dictionary_input": False,  # TODO-Dan Please document what each setting does
        "strictcheckvoltagesolution": False,
        "usevoltsolution": False,
        "fixcs": False,
        "skipjabr": False,
        "cutplane": False,
        "usemipstart": False,
        "useactivelossineqs": False,
        "useconvexformulation": False,
        "usemaxdispersion": False,
        "usemaxphasediff": False,
        "use_ef": False,
        "substitute_nonconv": False,
        "dopolar": False,
        "doslp_polar": False,
        "doac": False,
        "dodc": False,
        "doiv": False,
        "ivtype": None,
        "branchswitching_mip": False,
        "branchswitching_comp": False,
        "maxdispersion_deg": 0.0,
        "maxphasediff_deg": 360.0,
        "fixtolerance": 0.0,
    }

    if type(casefile) is not dict and casefile != "":
        settings["casefilename"] = casefile

    return settings


def get_default_graphics_settings(casefile):
    """Returns a dictionary holding all default settings"""

    settings = {
        "casefilename": None,
        "voltsfilename": None,  # TODO-Dan could you provide an example of how to use this?
        "gvfilename": None,  # TODO-Dan where is this used?
        "graphattrsfilename": None,  # TODO-Dan where is this used?
        "coordsfilename": None,
        "branchswitching_mip": False,
        "branchswitching_comp": False,
    }

    if type(casefile) is not dict and casefile != "":
        settings["casefilename"] = casefile

    return settings
