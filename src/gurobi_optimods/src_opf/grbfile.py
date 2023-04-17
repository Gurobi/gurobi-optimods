import sys
import os
import time
import math
import csv
import logging
import numpy as np

from .utils import (
    get_default_optimization_settings,
    get_default_graphics_settings,
    initialize_logger,
    remove_and_close_handlers,
    break_exit,
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
    alldata["usequadcostvar"] = False  # TODO-Dan What do we want to do with this?

    # TODO revisit default settings
    defaults = get_default_optimization_settings()
    read_settings_dict(alldata, settings, defaults)

    if alldata["doslp_polar"] and (int(alldata["dodc"]) + int(alldata["doiv"]) == 0):
        alldata["doac"] = True

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

    # TODO revisit default settings
    # DB this may be the wrong paradigm because we need the graphics
    # dict *before* alldata has been created, so as to pass it to the
    # function that creates alldata

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
        # Jarek, if we want to allow plotting via coordinates, we need to restructure this
    else:
        settings = get_default_optimization_settings()

    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    # Read all lines of the casefile
    logger.info("Building settings dictionary.")
    settings_dict = read_settings_build_dict(settings, lines)
    endtime = time.time()
    logger.info("Reading and building time: %f s." % (endtime - starttime))
    logger.info("")
    logger.info("Non-default settings read from file.")
    for s in settings_dict.keys():
        logger.info("  {} {}".format(s, settings_dict[s]))

    logger.info("")

    remove_and_close_handlers(logger, handlers)

    return settings_dict


def read_coordsettings_file(filename):
    coordsettings = {}
    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    linenum = 0
    # Read lines of configuration file and save options
    while linenum < len(lines):
        thisline = lines[linenum].split()

        if len(thisline) <= 0:
            linenum += 1
            continue
        if thisline[0][0] == "#":
            linenum += 1
            continue

        if thisline[0] == "END":
            break

        if thisline[0] == "coordsfilename":
            coordsettings["coordsfilename"] = thisline[1]
            print("COORDS at %s\n" % coordsettings["coordsfilename"])
            linenum += 1
            continue
        else:
            raise ValueError("Illegal option string %s." % thisline[0])

    # break_exit('poooooo')
    return coordsettings


def read_coords_build_dict(filename):
    """Reads csv file with bus coordinates, puts it into dict"""

    logger = logging.getLogger("OpfLogger")
    logger.info("Reading csv coordinates file %s and build dictionary." % filename)

    f = open(filename, "r")
    csvf = csv.reader(f)
    thelist = list(csvf)
    f.close()

    coords_dict = {}

    for line in range(1, len(thelist)):
        # print(thelist[line][0], float(thelist[line][3]), float(thelist[line][4]))
        coords_dict[line] = (float(thelist[line][3]), float(thelist[line][4]))

    print("Done reading coordinates\n")

    return coords_dict


def grbmap_coords_from_dict(alldata, coords_dict):
    """Maps with bus coordinates"""

    logger = logging.getLogger("OpfLogger")
    numbuses = alldata["numbuses"]
    buses = alldata["buses"]

    for bnum in range(1, numbuses + 1):
        # print(thelist[line][0], float(thelist[line][3]), float(thelist[line][4]))
        buses[bnum].lat = coords_dict[bnum][0]
        buses[bnum].lon = coords_dict[bnum][1]

    break_exit("mapped coords")


def read_settings_build_dict(settings, lines):
    """
    Helper function to read thru all lines of a previously given settings
    file and return a dictionary holding all non-default settings

    Parameters
    ----------
    settings : dictionary
        Settings dictionary provided by the user
    lines : list
        List of all lines of a previously read in settings file

    Returns
    -------
    dictionary
        A dictionary holding all non-default settings
    """

    settings_dict = {}
    linenum = 0
    # Read lines of configuration file and save options
    while linenum < len(lines):
        thisline = lines[linenum].split()

        if len(thisline) <= 0:
            linenum += 1
            continue

        if thisline[0][0] == "#":
            linenum += 1
            continue

        if thisline[0] == "END":
            break

        if thisline[0] not in settings.keys():
            raise ValueError("Illegal option string %s." % thisline[0])

        if len(thisline) > 2 and thisline[0] not in [
            "voltsfilename",
            "usevoltsolution",
        ]:
            raise ValueError("Illegal option string %s." % thisline[3])

        if len(thisline) == 2:
            settings_dict[thisline[0]] = thisline[1]

        elif len(thisline) == 1:
            settings_dict[thisline[0]] = True

        # handle special settings
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

            if settings_dict["usevoltsolution"]:
                raise ValueError(
                    "Illegal option combination. Cannot use both options voltsfilename and voltsolution."
                )

            settings_dict["voltsfilename"] = thisline[1]
            settings_dict["fixtolerance"] = float(thisline[2])

        elif thisline[0] == "usevoltsolution":
            if settings_dict["voltsfilename"] != None:
                raise ValueError(
                    "Illegal option combination. Cannot use both voltsfilename and voltsolution."
                )

            settings_dict["usevoltsolution"] = True  # forcing solution in casefile

        elif thisline[0] == "useconvexformulation":
            settings_dict["useconvexformulation"] = True
            settings_dict["doac"] = True

        elif thisline[0] == "doiv":
            settings_dict["doiv"] = True
            settings_dict["use_ef"] = True  # needed
            if len(thisline) > 1:
                settings_dict["ivtype"] = thisline[1]

        elif thisline[0] == "dographics":
            settings_dict["dographics"] = True
            if len(thisline) > 1:
                settings_dict["graphattrsfilename"] = thisline[1]

        elif thisline[0] == "coordsfilename":
            settings_dict["coordsfilename"] = True
            if len(thisline) > 1:
                settings_dict["coordsfilename"] = thisline[1]

        linenum += 1

    return settings_dict


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

        if len(thisline) <= 0:
            linenum += 1
            continue

        if thisline[0] == "END":
            break

        thresh[numfeatures] = int(thisline[0])
        colorstring[numfeatures] = thisline[1]
        sizeval[numfeatures] = int(thisline[2])
        numfeatures += 1

        linenum += 1

    logger.info("Read %d graphical features." % numfeatures)
    alldata["graphical"]["numfeatures"] = numfeatures
    alldata["graphical"]["thresh"] = thresh
    alldata["graphical"]["colorstring"] = colorstring
    alldata["graphical"]["sizeval"] = sizeval


def grbread_coords(alldata):
    """
    Function to read a csv file with bus coordinates for plotting.
    The csv filename is given via the coordsfilename setting

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    """

    logger = logging.getLogger("OpfLogger")
    numbuses = alldata["numbuses"]
    buses = alldata["buses"]
    IDtoCountmap = alldata["IDtoCountmap"]

    filename = alldata["coordsfilename"]
    logger.info("Reading csv coordinates file %s." % filename)

    f = open(filename, "r")
    csvf = csv.reader(f)
    thelist = list(csvf)
    f.close()

    for line in range(1, len(thelist)):
        # print(thelist[line][0], float(thelist[line][3]), float(thelist[line][4]))
        buses[line].lat = float(thelist[line][3])
        buses[line].lon = float(thelist[line][4])


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
    numbuses = alldata["numbuses"]
    buses = alldata["buses"]
    IDtoCountmap = alldata["IDtoCountmap"]

    filename = alldata["voltsfilename"]
    logger.info("reading volts file %s." % filename)

    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    linenum = 0

    while linenum < len(lines):
        thisline = lines[linenum].split()
        if len(thisline) > 0:

            if thisline[0] == "END":
                break
            elif thisline[0] == "bus":
                busid = int(thisline[1])
                buscount = IDtoCountmap[busid]

                vm = float(thisline[3])
                vadeg = float(thisline[5])
                varad = vadeg * math.pi / 180.0

                buses[buscount].inputvoltage = 1
                buses[buscount].inputV = vm
                buses[buscount].inputA_rad = varad
            else:
                raise ValueError(
                    "Illegal line %d in voltsfile %s." % (linenum + 1, filename)
                )

        linenum += 1

    logger.info("Read volts file.\n")
