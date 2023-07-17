import csv
import logging
import math

from gurobi_optimods.opf.utils import (
    check_settings_for_correct_type,
    get_default_optimization_settings,
)

logger = logging.getLogger(__name__)


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
):
    settings = dict(
        dopolar=bool(polar),
        use_ef=bool(useef),
        skipjabr=not bool(usejabr),
        usemipstart=bool(usemipstart),
        minactivebranches=float(minactivebranches),
        useactivelossineqs=bool(useactivelossineqs),
    )

    opftype = opftype.lower()
    if opftype in ["ac", "dc", "iv"]:
        settings["do" + opftype] = True
    else:
        raise ValueError(f"Unknown OPF type {opftype}")

    ivtype = ivtype.lower()
    if ivtype in ["plain", "aggressive"]:
        settings["ivtype"] = ivtype
    else:
        raise ValueError(f"Unknown ivtype {ivtype}")

    if branchswitching == 0:
        settings["branchswitching_mip"] = False
        settings["branchswitching_comp"] = False
    elif branchswitching == 1:
        settings["branchswitching_mip"] = True
        settings["branchswitching_comp"] = False
    elif branchswitching == 2:
        settings["branchswitching_mip"] = False
        settings["branchswitching_comp"] = True
    else:
        raise ValueError(f"Unknown branchswitching setting {branchswitching}.")

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
