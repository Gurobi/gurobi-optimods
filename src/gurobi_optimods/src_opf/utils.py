from enum import Enum
import logging
import sys


def break_exit(foo):
    """
    For Dan only
    Will be definitely removed in the final version
    """

    stuff = input("(" + foo + ") break> ")
    if stuff == "x" or stuff == "q":
        sys.exit("bye")


def initialize_logger(loggername, logfile="", usefilehandler=False):
    logger = logging.getLogger("OpfLogger")
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(message)s")
    stdouthandler = logging.StreamHandler(stream=sys.stdout)
    stdouthandler.setLevel(logging.INFO)
    stdouthandler.setFormatter(formatter)
    handlers = [stdouthandler]
    if usefilehandler and logfile:
        filehandler = logging.FileHandler(filename=logfile)
        filehandler.setLevel(logging.INFO)
        filehandler.setFormatter(formatter)
        handlers.append(filehandler)
    for h in handlers:
        logger.addHandler(h)

    return logger, handlers


def remove_and_close_handlers(logger, handlers):
    for h in handlers:
        logger.removeHandler(h)
        h.close()


class OpfType(Enum):
    AC = "AC"
    DC = "DC"
    IV = "IV"


def get_default_optimization_settings():
    """Returns a dictionary holding all default settings"""

    settings = {
        "casefilename": None,
        "voltsfilename": None,  # TODO-Dan could you provide an example of how to use this? I will.  It is a plain text file that has, for each bus, a line of the form "bus 8 M 1.099999e+00 A 9.051238e-01" (bus numbers 1 through N) plus a final "END" line
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
        "voltsfilename": None,  # TODO-Dan could you provide an example of how to use this? Will do
        "gvfilename": None,  # TODO-Dan where is this used?  It may only get used if we are using sfdp
        "graphattrsfilename": None,  # TODO-Dan where is this used? These are color choices for rendering graphics
        "coordsfilename": None,
        "branchswitching_mip": False,
        "branchswitching_comp": False,
    }

    return settings
