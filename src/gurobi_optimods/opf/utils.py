import logging
import sys
from enum import Enum


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
        "doac": False,
        "dodc": False,
        "doiv": False,
        "dopolar": False,
        # "doslp_polar": False, # Not used yet, this is for future work
        "use_ef": True,
        "skipjabr": False,
        "ivtype": "aggressive",
        "branchswitching_mip": False,
        # Formulation for branch-switching where the binary variables simply multiply the continuous variables.
        # Sometimes it works better than branchswitching_mip. Only applicable for AC
        "branchswitching_comp": False,
        "usemipstart": False,
        "minactivebranches": 0.95,
        # New linear inequalities developed and implemented by Dan.
        # They are outer approximations of the JABR inequalities
        "useactivelossineqs": False,
        #############################
        # The following settings should currently not be disclosed
        # For now keep for us, mainly used for debugging and experimenting with heuristics
        "fixcs": False,  # (approximately) fix c, s variables if a voltage solution was read in
        "fixtolerance": 1.0e-5,
        # Heuristics to help NL solver find a good solution
        "usemaxdispersion": False,  # difference between all bus angles is small
        "usemaxphasediff": False,  # difference between 2 adjacent branches is small
        "maxdispersion_deg": 0.0,
        "maxphasediff_deg": 360.0,
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
                "lpfilename",
                "gurobiparamfile",
                "ivtype",
            ]
            and settingval != None
            and not isinstance(settingval, str)
        ):
            wrongsetting = settingname
            settingtype = "String"
            break
        # Check boolean settings
        if settingname in [
            "fixcs",
            "skipjabr",
            "usemipstart",
            "useactivelossineqs",
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
