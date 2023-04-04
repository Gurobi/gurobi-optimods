from enum import Enum
import sys


def break_exit(foo):
    """
    For Dan only
    Will be definitely removed in the final version
    """

    stuff = input("(" + foo + ") break> ")
    if stuff == "x" or stuff == "q":
        sys.exit("bye")


class OpfType(Enum):
    AC = "AC"
    DC = "DC"
    IV = "IV"


def get_default_optimization_settings():
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

    return settings


def get_default_graphics_settings():
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

    return settings
