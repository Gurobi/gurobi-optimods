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


    Parameters
    ----------
    logfile : string, optional
        Name of log file. Can be empty

    Returns
    -------
    dictionary
        Returns a dictionary with a few initialized default fields
    """
    alldata = {}
    alldata["LP"] = {}
    alldata["MIP"] = {}
    alldata["logfile"] = logfile

    return alldata


def read_optimization_settings(alldata, settings):
    """
    Reads settings dictionary holding settings for an optimization call
    and saves non-default settings to alldata dictionary.
    Also performs basic checks for incompatible settings


    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    settings : dictionary
        Dictionary holding settings used for an optimization call provided by the user
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
        logger.info("  {} {}".format(s, alldata[s]))

    logger.info("")


def read_graphics_settings(alldata, settings):
    """
    Reads settings dictionary holding settings for a graphics call
    and saves non-default settings to alldata dictionary

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    settings : dictionary
        Dictionary holding settings used for a graphics call provided by the user
    """

    defaults = get_default_graphics_settings()
    read_settings_dict(alldata, settings, defaults)


def read_settings_dict(alldata, inputsettings, defaults):
    """
    Sets input settings in alldata from a settings dict.

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    inputsettings : dictionary
        Dictionary holding settings provided by the user
    defaults: dictionary
        Dictionary holding optimization or graphics default settings
    """

    for s in inputsettings.keys():
        if s not in defaults.keys():
            raise ValueError("Illegal option string %s." % s)
        # Possible TODO do we have to check for correct types of values inputsettings[s]?
        defaults[s] = inputsettings[s]

    check_settings_for_correct_type(defaults)

    for s in defaults.keys():
        alldata[s] = defaults[s]


def read_settings_file(filename, graphics=False):
    """
    Function to read settings from a file

    Parameters
    ----------
    filename : string
        Path to text based file holding user settings
    graphics : boolean, optional
        If set to true, then graphic specific settings are expected. Otherwise,
        settings for an optimization call are expected

    Returns
    -------
    dictionary
        A dictionary holding all non-default settings
    """

    logger, handlers = initialize_logger("SettingsReadingLogger")

    logger.info("Reading settings file %s." % filename)
    starttime = time.time()
    settings = {}
    if graphics:
        settings = get_default_graphics_settings()
    else:
        settings = get_default_optimization_settings()

    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    # Read all lines of the casefile
    logger.info("Building settings dictionary.")
    settings_dict = read_settings_build_dict(settings, lines)
    endtime = time.time()
    logger.info(
        "Settings reading and dictionary building time: %f s." % (endtime - starttime)
    )
    logger.info("")
    logger.info("User set settings read from file.")
    for s in settings_dict.keys():
        logger.info("  {} {}".format(s, settings_dict[s]))

    logger.info("")

    remove_and_close_handlers(logger, handlers)

    return settings_dict


def read_settings_build_dict(settings, lines):
    """
    Helper function to read thru all lines of a previously given settings
    file and return a dictionary holding all non-default settings

    Parameters
    ----------
    settings : dictionary
        Settings dictionary currently holding default settings which will
        be overwritten by settings defined by lines
    lines : list
        List of all lines of a previously read in settings file

    Returns
    -------
    dictionary
        A dictionary holding all user set settings
    """

    settings_dict = {}
    # Read lines of configuration file and save options
    for line in lines:
        thisline = line.split()
        # Skip empty lines and comments
        if len(thisline) <= 0 or thisline[0][0] == "#":
            continue
        # End of settings
        if thisline[0] == "END":
            break
        # Check whether we know the setting
        if thisline[0] not in settings.keys():
            raise ValueError("Illegal option string %s." % thisline[0])
        # We only have one setting which accepts 3 inputs in one line
        if len(thisline) > 2 and thisline[0] not in [
            "voltsfilename",
        ]:
            raise ValueError("Illegal option string %s." % thisline[3])
        # Set setting with 2 inputs
        if len(thisline) == 2:
            if thisline[0] in ["maxdispersion_deg", "maxphasediff_deg", "fixtolerance"]:
                settings_dict[thisline[0]] = float(thisline[1])
            else:
                settings_dict[thisline[0]] = thisline[1]
        # Setting with 1 input are booleans, if they are present then they are True
        elif len(thisline) == 1:
            settings_dict[thisline[0]] = True
        # Handle special settings which may or may not take a second argument
        if thisline[0] == "usemaxdispersion":
            settings_dict["usemaxdispersion"] = True
            if len(thisline) > 1:
                settings_dict["maxdispersion_deg"] = float(thisline[1])

        elif thisline[0] == "usemaxphasediff":
            settings_dict["usemaxphasediff"] = True
            if len(thisline) > 1:
                settings_dict["maxphasediff_deg"] = float(thisline[1])

        elif thisline[0] == "strictcheckvoltagesolution":
            settings_dict["strictcheckvoltagesolution"] = True
            if len(thisline) > 1:
                settings_dict["voltsfilename"] = thisline[1]

        elif thisline[0] == "voltsfilename":
            if len(thisline) < 3:
                raise ValueError(
                    "Voltsfilename option requires filename and tolerance."
                )

            settings_dict["voltsfilename"] = thisline[1]
            settings_dict["fixtolerance"] = float(thisline[2])

        elif thisline[0] == "useconvexformulation":
            settings_dict["useconvexformulation"] = True
            settings_dict["doac"] = True

        elif thisline[0] == "doiv":
            settings_dict["doiv"] = True
            settings_dict["use_ef"] = True  # needed
            if len(thisline) > 1:
                settings_dict["ivtype"] = thisline[1]

    return settings_dict


def read_coords_file_csv(filename):
    """
    Function to read bus coordinates from a csv file

    Parameters
    ----------
    filename : string
        Path to text based file holding bus coordinates

    Returns
    -------
    dictionary
        A dictionary holding the respective coordinates
    """

    logger, handlers = initialize_logger("CoordsReadingLogger")
    logger.info("Reading csv coordinates file %s and building dictionary." % filename)

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

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    coords_dict : dictionary
        Dictionary holding a mapping of bus index to a given coordinate
    """

    numbuses = alldata["numbuses"]
    buses = alldata["buses"]

    for bnum in range(1, numbuses + 1):
        index = buses[bnum].nodeID
        buses[bnum].lat = coords_dict[index][0]
        buses[bnum].lon = coords_dict[index][1]


def grbread_graphattrs(alldata, filename):
    """
    Function to read graphical attributes from a user-given file and
    save these into the alldata dictionary

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    filename : string
        Path to text based file holding user-given graph attributes
    """

    logger = logging.getLogger("OpfLogger")
    logger.info("Reading graphical attributes file %s." % filename)
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

    logger.info("Succesfully read %d graphical features." % numfeatures)
    alldata["graphical"]["numfeatures"] = numfeatures
    alldata["graphical"]["thresh"] = thresh
    alldata["graphical"]["colorstring"] = colorstring
    alldata["graphical"]["sizeval"] = sizeval


def grbreadvoltsfile(alldata):
    """
    TODO-Dan Review description
    Function to read a file holding input voltages for buses.
    The voltage filename is given via the voltsfilename setting

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    """

    logger = logging.getLogger("OpfLogger")
    buses = alldata["buses"]
    IDtoCountmap = alldata["IDtoCountmap"]

    filename = alldata["voltsfilename"]
    logger.info("Reading volts file %s." % filename)

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
            raise ValueError(
                "Illegal line %d in voltsfile %s." % (linenum + 1, filename)
            )

        linenum += 1

    logger.info("Succesfully read volts file.\n")
