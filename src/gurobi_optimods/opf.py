import logging

import gurobipy as gp

from .src_opf.grbcasereader import (
    read_case,
    read_case_file,
    read_case_file_mat,
    turn_opf_dict_into_mat_file,
)

from .src_opf.grbfile import (
    initialize_data_dict,
    read_settings_file,
    read_optimization_settings,
    read_graphics_settings,
    read_coords_file_csv,
    grbmap_coords_from_dict,
    grbread_graphattrs,
)
from .src_opf.grbformulator import construct_and_solve_model
from .src_opf.grbgraphical import generate_solution_figure
from .src_opf.utils import initialize_logger, remove_and_close_handlers


def solve_opf_model(settings, case, logfile=""):
    """
    Constructs an OPF model from given data and solves it with Gurobi.
    Returns a result dictionary following MATPOWER notation

    :param settings: Dictionary holding settings
    :type settings: dictionary
    :param case: Dictionary holding case data
    :type case: dictionary
    :param logfile: Name of log file, defaults to ""
    :type logfile: str, optional

    :return: A result dictionary following MATPOWER notation, e.g., the objective value
             is at solution["f"] if solution["success"] = 1
    :rtype: dictionary
    """

    # Initialize output and file handler and start logging
    logger, handlers = initialize_logger("OpfLogger", logfile, True)

    # Initilize data dictionary
    alldata = initialize_data_dict(logfile)

    # Read settings file/dict and save them into the alldata dict
    read_optimization_settings(alldata, settings)

    # Read case file/dict and populate the alldata dictionary
    read_case(alldata, case)

    # Construct and solve model using given case data and user settings
    solution = construct_and_solve_model(alldata)

    # Remove and close all logging handlers
    remove_and_close_handlers(logger, handlers)

    return solution


def generate_opf_solution_figure(settings, case, coords, solution):
    """
    Reads the given case and returns a plotly figure object.
    Ideally the solution has been computed by the `solve_opf_model` function

    :param settings: Dictionary holding user settings
    :type settings: dictionary
    :param case: Dictionary holding case data
    :type case: dictionary
    :param coords: Dictionary holding bus coordinates
    :type coords: dictionary
    :param solution: Dictionary holding solution data following the MATPOWER notation as returned
                      by the solve_opf_model function
    :type solution: dictionary


    :return: A plotly figure objects which can be displayed via the show() function,
             see https://plotly.com/python-api-reference/generated/plotly.graph_objects.Figure.html
    :rtype: :class:`plotly.graph_objects.Figure`
    """

    # Initialize output and file handler and start logging
    logger, handlers = initialize_logger("OpfLogger")

    # Initilize data dictionary
    alldata = initialize_data_dict()

    # Read settings file/dict and save them into the alldata dict
    read_graphics_settings(alldata, settings)

    # Read case file/dict and populate the alldata dictionary
    read_case(alldata, case)

    # Special settings for graphics
    alldata["graphical"] = {}
    alldata["graphical"]["numfeatures"] = 0
    if alldata["graphattrsfilename"] != None:
        grbread_graphattrs(alldata, alldata["graphattrsfilename"])

    # Map given coordinate data to network data
    grbmap_coords_from_dict(alldata, coords)

    # Generate a plotly figure object representing the given solution for the network
    fig = generate_solution_figure(alldata, solution)

    # Remove and close all logging handlers
    remove_and_close_handlers(logger, handlers)

    return fig


def read_settings_from_file(settingsfile, graphics=False):
    """
    Helper function for users.
    Constructs a settings dictionary which can be used as input
    for all other OPF API functions

    :param settingsfile: Name of and possibly full path to settings file
    :type settingsfile: string
    :param graphics: `True` is the settingsfile holds only graphics settings,
                     `False` otherwise
    :type graphics: bool, optional

    :return: Dictionary object of user settings read in from settingsfile
              which can be used in other OPF API functions
    :rtype: dictionary
    """

    settings_dict = read_settings_file(settingsfile, graphics)

    return settings_dict


def read_case_from_file(casefile):
    """
    Helper function for users.
    Reads a `.m` text file and constructs a case dictionary,
    which can be used as input for other OPF API functions

    :param casefile: Name of and possibly full path to case file given as `.m` file
                     The `.m` file should be in standard MATPOWER notation
    :type casefile: str

    :return: Dictionary object of the given case which can be used in
             other OPF API functions
    :rtype: dictionary
    """

    case_dict = read_case_file(casefile)

    return case_dict


def read_case_from_mat_file(casefile):
    """
    Helper function for users.
    Reads a `.mat` data file and constructs a case dictionary,
    which can be used as input for other OPF API functions

    :param casefile: Name of and possibly full path to case file given as `.mat` file
                     The `.mat` file should be in standard MATPOWER notation
    :type casefile: str

    :return: Dictionary object of the given case which can be used in
             other OPF API functions
    :rtype: dictionary
    """

    case_dict = read_case_file_mat(casefile)

    return case_dict


def turn_solution_into_mat_file(solution, matfilename=""):
    """
    Writes a `.mat` data file in MATPOWER notation out of an OPF solution dictionary

    :param solution: OPF solution dictionary
    :type solution: dictionary
    :param matfilename: Name of `.mat` file where to write the solution data,
                        defaults to `result.mat`
    :type matfilename: str, optional
    """

    # Initialize output and file handler and start logging
    logger, handlers = initialize_logger("OpfLogger")

    # Set a default output file name
    if matfilename == "":
        matfilename = "result.mat"

    # Check for .mat suffix
    if not matfilename[-4:] == ".mat":
        matfilename += ".mat"
    logger.info(
        f"Generating .mat file {matfilename} out of given OPF solution dictionary."
    )

    # Generate .mat file out of solution
    turn_opf_dict_into_mat_file(solution, matfilename)

    # Remove and close all logging handlers
    remove_and_close_handlers(logger, handlers)


def read_coords_from_csv_file(coordsfile):
    """
    Helper function for users.
    Construct a coordinate dictionary which can be used as input
    for the `generate_opf_solution_figure` function

    :param coordsfile: Name of and possibly full path to case file given as `.csv` file
    :type coordsfile: str

    :return: Dictionary of the given coordinates which can be used in
             the `generate_opf_solution_figure` function
    :rtype: dictionary
    """

    coord_dict = read_coords_file_csv(coordsfile)

    return coord_dict
