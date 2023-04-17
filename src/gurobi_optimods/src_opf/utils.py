from enum import Enum
import logging
import sys  # TODO remove in final version


def break_exit(foo):
    """
    For Dan only
    Will be definitely removed in the final version
    """

    stuff = input("(" + foo + ") break> ")
    if stuff == "x" or stuff == "q":
        sys.exit("bye")


def initialize_logger(loggername, logfile="", usefilehandler=False):
    """
    Initializes a logger object from the logging package


    Parameters
    ----------
    loggername : string
        Name of logger

    logfile : string, optional
        Name of log file. Can be empty

    usefilehandler: boolean, optional
        If set to True, then an aditional filehandler is created and the
        returned logger will write to a file specified by logfile


    Returns
    -------
    logging.logger, list
        Returns the created logger object and a list of associated handler,
        a StreamHandler and possibly a FileHandler
    """
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
    """
    Removes and closes all handler associated to a given logger object


    Parameters
    ----------
    logger : logging.logger
        Logger object

    handlers : list
        List of handlers associated to the logger object
    """

    for h in handlers:
        logger.removeHandler(h)
        h.close()


class OpfType(Enum):
    """
    Defines possible OPF formulation types
    """

    AC = "AC"
    DC = "DC"
    IV = "IV"


def get_default_optimization_settings():
    """
    Returns a dictionary holding default settings for an optimization call

    Returns
    -------
    dictionary
        A dictionary holding default setting for an optimization call
    """

    settings = {
        "coords_dict": None,
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
    """
    Returns a dictionary holding default settings for a graphics call

    Returns
    -------
    dictionary
        A dictionary holding default setting for an graphics call
    """

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
