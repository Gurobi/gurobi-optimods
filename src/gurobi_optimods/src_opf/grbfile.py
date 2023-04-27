import time
import math
import csv
import logging

from .utils import (
    get_default_optimization_settings,
    get_default_graphics_settings,
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
        logger.info(f"  {s} {alldata[s]}")

    logger.info("")


def read_graphics_settings(alldata, settings):
    """
    Reads settings dictionary for a graphics call
    and saves all settings to alldata dictionary

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param settings: Dictionary holding settings used for a graphics call provided by the user
    :type settings: dict
    """

    defaults = get_default_graphics_settings()
    read_settings_dict(alldata, settings, defaults)


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


def read_coords_file_csv(filename):
    """
    Reads bus coordinates from a `.csv` file

    :param filename: Path to text based file holding bus coordinates
    :type filename: string

    :return: Dictionary holding the respective coordinates
    :rtype: dict
    # TODO-Dan what are the other inputs of the csv file?
    """

    logger, handlers = initialize_logger("CoordsReadingLogger")
    logger.info(f"Reading csv coordinates file {filename} and building dictionary.")

    with open(filename, mode="r") as infile:
        reader = csv.reader(infile)
        coords_dict = {
            int(rows[1]): (
                rows[3],
                rows[4],
            )
            for rows in reader
            if rows[0] != "index"
        }

    logger.info("Done reading coordinates.\n")
    remove_and_close_handlers(logger, handlers)

    return coords_dict


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


def grbread_graphattrs(alldata, filename):
    """
    Reads graphical attributes from a user-given file and
    saves these into the alldata dictionary

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param filename: Path to text based file holding user-given graph attributes
    :type filename: str
    """

    logger = logging.getLogger("OpfLogger")
    logger.info(f"Reading graphical attributes file {filename}.")
    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    numfeatures = 0

    thresh = {}
    colorstring = {}
    sizeval = {}

    linenum = 0
    # Read lines of configuration file and save options
    while linenum < len(lines):
        thisline = lines[linenum].split()
        # Skip empty lines and comments
        if len(thisline) <= 0 or thisline[0] == "#":
            linenum += 1
            continue
        # End of graph attributes
        if thisline[0] == "END":
            break

        thresh[numfeatures] = int(thisline[0])
        colorstring[numfeatures] = thisline[1]
        sizeval[numfeatures] = int(thisline[2])
        numfeatures += 1

        linenum += 1

    logger.info(f"Succesfully read {numfeatures} graphical features.")
    alldata["graphical"]["numfeatures"] = numfeatures
    alldata["graphical"]["thresh"] = thresh
    alldata["graphical"]["colorstring"] = colorstring
    alldata["graphical"]["sizeval"] = sizeval


def grbreadvoltsfile(alldata):
    """
    TODO-Dan Review description
    Read a file holding input voltages for buses and saves
    necessary data into the alldata dictionary

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict

    :raises ValueError: Illegal line in voltsfile
    """

    logger = logging.getLogger("OpfLogger")
    buses = alldata["buses"]
    IDtoCountmap = alldata["IDtoCountmap"]

    filename = alldata["voltsfilename"]
    logger.info(f"Reading volts file {filename}.")

    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    linenum = 0
    while linenum < len(lines):
        thisline = lines[linenum].split()
        # Skip empty lines and comments
        if len(thisline) <= 0 or thisline[0] == "#":
            linenum += 1
            continue
        # End of voltage information
        if thisline[0] == "END":
            break

        if thisline[0] == "bus":
            busid = int(thisline[1])
            buscount = IDtoCountmap[busid]

            vadeg = float(thisline[5])
            buses[buscount].inputvoltage = 1
            buses[buscount].inputV = float(thisline[3])  # Vm
            buses[buscount].inputA_rad = vadeg * math.pi / 180.0  # VaRad
        else:
            raise ValueError(f"Illegal line {linenum + 1} in voltsfile {filename}.")

        linenum += 1

    logger.info("Succesfully read volts file.\n")
