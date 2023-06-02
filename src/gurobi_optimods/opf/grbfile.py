import time
import math
import csv
import logging

from gurobi_optimods.opf.utils import (
    get_default_optimization_settings,
    initialize_logger,
    remove_and_close_handlers,
    check_settings_for_correct_type,
)


def initialize_data_dict(logfile=""):
    """
    Initializes a dictionary holding all necessary data

    :param logfile: Name of log file, defaults to ""
    :type logfile: str, optional

    :return: Returns a dictionary with a few initialized default fields
    :rtype: dict
    """
    alldata = {}
    alldata["LP"] = {}  # continuous variables
    alldata["MIP"] = {}  # discrete variables
    alldata["logfile"] = logfile

    return alldata


def construct_settings_dict(
    opftype,
    polar,
    useef,
    usejabr,
    ivtype,
    branchswitching,
    usemipstart,
    minactivebranches,
    useactivelossineqs,
    additional_settings,
):
    """
    Initializes a dictionary holding all necessary data

    :param logfile: Name of log file, defaults to ""
    :type logfile: str, optional
    :param opftype: String telling the desired OPF model type. Available are `AC`, `DC`, `IV`, defaults to `AC`
    :type opftype: str, optional
    :param polar: Controls whether polar formulation should be used, defaults to `False`
    :type polar: bool, optional
    :param useef: Controls whether bilinear variables e, f and corresponding constraints should be used, defaults to `True`.
                  Has only an effect if ``opftype`` equals `AC`
    :type useef: bool, optional
    :param usejabr: Controls whether JABR inequalities should be added, defaults to `True`.
                    Has only an effect if ``opftype`` equals `AC`
    :type usejabr: bool, optional
    :param ivtype: States what type of IV formulation should be used. Availale are `aggressive` and `plain`,
                   defaults to `aggressive`
    :type ivtype: str, optional
    :param branchswitching: Controls whether discrete variable for turning on/off branches should be used.
                            0 = don't use discrete variables (all branches are on)
                            1 = use binary variables to constrain bounds of (re)active power injections
                            2 = use binary variables and multiply them with the (re)active power injections
                            Usually, setting 1 works better than 2. Defaults to 0
    :type branchswitching: int, optional
    :param usemipstart: Controls whether a pre-defined MIPStart should be used. Has only an effect if
                        branchswitching > 0. Deftault to `True`
    :type usemipstart: bool, optional
    :param additional_settings: Dictionary holding additional settings. For more details, please refer to
                                the documentation. Defaults to an empty dictionary
    :type additional_settings: dict, optional

    :return: Returns a dictionary with a few initialized default fields
    :rtype: dict
    """

    settings = {}

    opftypeL = opftype.lower()
    if opftypeL not in ["ac", "dc", "iv"]:
        raise ValueError(f"Unknown OPF type {opftype}.")

    settings["do" + opftypeL] = True
    settings["dopolar"] = polar
    settings["use_ef"] = useef
    settings["skipjabr"] = not usejabr

    ivtypeL = ivtype.lower()
    if ivtypeL not in ["plain", "aggressive"]:
        raise ValueError(f"Unknown ivtype {ivtype}.")

    settings["ivtype"] = ivtypeL
    if branchswitching == 0:
        settings["branchswitching_mip"] = settings["branchswitching_comp"] = False
    elif branchswitching == 1:
        settings["branchswitching_mip"] = True
        settings["branchswitching_comp"] = False
    elif branchswitching == 2:
        settings["branchswitching_mip"] = False
        settings["branchswitching_comp"] = True
    else:
        raise ValueError(f"Unknown branchswitching setting {branchswitching}.")

    settings["usemipstart"] = usemipstart
    settings["minactivebranches"] = minactivebranches

    settings["useactivelossineqs"] = useactivelossineqs

    for s in additional_settings:
        settings[s] = additional_settings[s]

    return settings


def read_optimization_settings(alldata, settings):
    """
    Reads settings dictionary for an optimization call
    and saves all settings to alldata dictionary.

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param settings: Dictionary holding settings used for an optimization call provided by the user
    :type settings: dict

    :raise ValueError: Incompatible settings
    """

    # Currently Hard-coded
    # We leave it like this for now, it needs further experimentation in a future release
    alldata["usequadcostvar"] = False

    defaults = get_default_optimization_settings()
    read_settings_dict(alldata, settings, defaults)

    if int(alldata["doac"]) + int(alldata["dodc"]) + int(alldata["doiv"]) != 1:
        raise ValueError(
            "Illegal option combination. Have to use exactly 1 of options [doac, dodc, doiv]."
        )

    logger = logging.getLogger("OpfLogger")
    logger.info("All settings:")
    for s in defaults.keys():
        if s == "skipjabr":
            logger.info(f"  usejabr {not alldata[s]}")
        # Don't print hidden settings
        elif s not in [
            "fixcs",
            "fixtolerance",
            "usemaxdispersion",
            "usemaxphasediff",
            "maxdispersion_deg",
            "maxphasediff_deg",
        ]:
            logger.info(f"  {s} {alldata[s]}")

    logger.info("")


def read_settings_dict(alldata, inputsettings, defaults):
    """
    Sets input settings in alldata from a settings dict.
    Also checks for correct data types of settings

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param inputsettings: Dictionary holding settings provided by the user
    :type inputsettings: dict
    :param defaults: Dictionary holding optimization or graphics default settings
    :type defaults: dict

    :raises ValueError: Unknown option string
    """

    for s in inputsettings.keys():
        if s not in defaults.keys():
            raise ValueError(f"Unknown option string {s}.")
        defaults[s] = inputsettings[s]

    check_settings_for_correct_type(defaults)

    for s in defaults.keys():
        alldata[s] = defaults[s]


def read_file_csv(filename, data):
    """
    Reads bus data from a `.csv` file

    :param filename: Path to text based file holding bus data
    :type filename: str
    :param data: Name of data to be read in
    :type data: str

    :return: Dictionary holding the respective data of the form busID : (data1, data2)
    :rtype: dict
    """

    logger, handlers = initialize_logger("CoordsReadingLogger")
    logger.info(f"Reading csv {data} file {filename} and building dictionary.")

    with open(filename, mode="r") as infile:
        reader = csv.reader(infile)
        # Sanity check
        data_dict = {}
        first = True
        for rows in reader:
            # Skip first row of the .csv file
            if first:
                first = False
                continue

            if len(rows) != 5:
                raise ValueError(
                    f"Incorrect input in {data} .csv file {filename}. Number of columns does not equal 5."
                )

            data_dict[int(rows[1])] = (
                float(rows[3]),
                float(rows[4]),
            )

    logger.info(f"Done reading {data}.\n")
    remove_and_close_handlers(logger, handlers)

    return data_dict


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

    buses = alldata["buses"]
    IDtoCountmap = alldata["IDtoCountmap"]

    for v in volts_dict:
        buscount = IDtoCountmap[v]
        buses[buscount].inputvoltage = 1
        buses[buscount].inputV = volts_dict[v][0]  # Vm
        vadeg = volts_dict[v][1]
        buses[buscount].inputA_rad = vadeg * math.pi / 180.0  # VaRad
