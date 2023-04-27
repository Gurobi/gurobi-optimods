import logging
import sys
from enum import Enum


def break_exit(foo):
    """
    For Dan only
    Will be definitely removed in the final version
    # TODO remove it
    """

    stuff = input("(" + foo + ") break> ")
    if stuff == "x" or stuff == "q":
        sys.exit("bye")


def initialize_logger(loggername, logfile="", usefilehandler=False):
    """
    Initializes a logger object from the logging package

    :param loggername: Name of logger
    :type loggername: str
    :param logfile: Name of log file, defaults to ""
    :type logfile: str, optional
    :param usefilehandler: If set to True, then an aditional filehandler is created and the
                           returned logger will write to a file specified by logfile, defaults
                           to False
    :type usefilehandler: bool, optional

    :rtype: :class: `logging.logger`, list
    :return: Returns the created logger object and a list of associated handlers,
             a :class: `logging.StreamHandler` and possibly a :class: `logging.FileHandler`
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

    :param logger: Logger object
    :type logger: :class: `logging.logger`
    :param handler: List of handlers associated to the logger object
    :type handlers: list
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

    :rtype: dict
    :return: Dictionary holding default setting for an optimization call
    """

    settings = {
        "lpfilename": None,
        "gurobiparamfile": None,
        "strictcheckvoltagesolution": False,  # will be removed
        "doac": False,  # TODO combine to a opftype argument
        "dodc": False,
        "doiv": False,
        "dopolar": False,
        # "doslp_polar": False, # Not used yet, this is for future work
        "use_ef": True,
        "skipjabr": False,
        "ivtype": "aggressive",
        "branchswitching_mip": False,
        # Formulation for branch-switching where the binary variables simply multiply the continuous variables.
        # Sometimes it works better. Only applicable for AC.
        "branchswitching_comp": False,
        "usemipstart": False,
        "useactivelossineqs": False,  # new linear inequalities developed and implemented by Dan
        # the following settings should currently not be disclosed
        # for now keep for us, mainly used for debugging and experimenting with heuristics
        "voltsfilename": None,  # TODO-Dan could you provide an example of how to use this? I will.  It is a plain text file that has, for each bus, a line of the form "bus 8 M 1.099999e+00 A 9.051238e-01" (bus numbers 1 through N) plus a final "END" line
        "fixcs": False,  # (approximately) fix c, s variables if a voltage solution was read in
        "fixtolerance": 1.0e-5,
        # heuristics to help NL solver find a good solution
        "usemaxdispersion": False,  # difference between all bus angles is small
        "usemaxphasediff": False,  # difference between 2 adjacent branches is small
        "maxdispersion_deg": 0.0,
        "maxphasediff_deg": 360.0,
    }

    return settings


def get_default_graphics_settings():
    """
    Returns a dictionary holding default settings for a graphics call

    :rtype: dict
    :return: Dictionary holding default setting for a graphics call
    """

    settings = {
        "voltsfilename": None,  # TODO-Dan could you provide an example of how to use this? Will do
        "graphattrsfilename": None,  # TODO-Dan where is this used? These are color choices for rendering graphics
    }

    return settings


def check_settings_for_correct_type(settings):
    """
    Checks whether the provided settings have correct types

    :param settings: Settings dictionary to be checked for correct data types
    :type settings: dict

    :raises ValueError: Wrong setting type
    """

    wrongsetting = None
    settingtype = ""
    for settingname in settings:
        settingval = settings[settingname]
        # Check string settings
        if (
            settingname
            in [
                "voltsfilename",
                "lpfilename",
                "gurobiparamfile",
                "ivtype" "graphattrsfilename",
            ]
            and settingval != None
            and not isinstance(settingval, str)
        ):
            wrongsetting = settingname
            settingtype = "String"
            break
        # Check boolean settings
        if settingname in [
            "strictcheckvoltagesolution",
            "fixcs",
            "skipjabr",
            "usemipstart",
            "useactivelossineqs",
            "useconvexformulation",
            "usemaxdispersion",
            "usemaxphasediff",
            "use_ef",
            "dopolar",
            "doac",
            "dodc",
            "doiv",
            "branchswitching_mip",
            "branchswitching_comp",
        ] and not isinstance(settingval, bool):
            # Accept 1 or 0
            if (isinstance(settingval, int) or isinstance(settingval, float)) and (
                settingval == 1 or settingval == 0
            ):
                settings[settingname] = True if settingval == 1 else False
            else:
                wrongsetting = settingname
                settingtype = "Boolean"
                break
        # Check numeric settings
        if (
            settingname in ["maxdispersion_deg", "maxphasediff_deg", "fixtolerance"]
            and not isinstance(settingval, int)
            and not isinstance(settingval, float)
        ):
            wrongsetting = settingname
            settingtype = "Int or Float"
            break

    if wrongsetting != None:
        raise ValueError(f"Setting {wrongsetting} is not of type {settingtype}.")
