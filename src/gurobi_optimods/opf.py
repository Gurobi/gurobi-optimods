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
    Construct an OPF model from given data and solve it with Gurobi


    Parameters
    ----------
    settings : dictionary
        Dictionary holding settings
    case : dictionary
        Dictionary holding case data
    logfile: string
        Name of log file. Can be empty


    Returns
    -------
    OrderedDict, float
        A feasible solution point if any was found as an OrderedDict and
        the final objective value
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
    solution, objval = construct_and_solve_model(alldata)

    # Remove and close all logging handlers
    remove_and_close_handlers(logger, handlers)

    return solution, objval


def generate_opf_solution_figure(settings, case, coords, solution, objval):
    """
    Read the given case and plot a given solution.
    In best case the solution has been computed by the solve_opf_model function
    so the indices of variables fit the ones computed by the previously
    constructed model.


    Parameters
    ----------
    settings : dictionary
        Dictionary holding settings
    case : dictionary
        Dictionary holding case data
    coords : dictionary
        Dictionary holding bus coordinates
    solution : OrderedDict
        OrderedDictionary holding solution data
    objval : float
        Objective value

    Returns
    -------
    plotly.graph_objects.Figure
        A plotly figure objects which can be displayed via the show() function,
        see https://plotly.com/python-api-reference/generated/plotly.graph_objects.Figure.html
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
    fig = generate_solution_figure(alldata, solution, objval)

    # Remove and close all logging handlers
    remove_and_close_handlers(logger, handlers)

    return fig


def read_settings_from_file(settingsfile, graphics=False):
    """
    Helper function for users
    Used to construct a settings dictionary which can be used as input
    for all other API functions

    Parameters
    ----------
    settingsfile : string
        Name of and possibly full path to settings file
    graphics : boolean, optional
        If set to true, then the function expects a settings file
        with special graphics settings


    Returns
    -------
    dictionary
        Dictionary object of the given settings which is to be used in
        other API functions
    """

    settings_dict = read_settings_file(settingsfile, graphics)

    return settings_dict


def read_case_from_file(casefile):
    """
    Helper function for users
    Used to construct a case dictionary which can be used as input
    for other API functions

    Parameters
    ----------
    casefile :
        Name of and possibly full path to case file given as .m file
        The .m file should be in standard MATPOWER notation

    Returns
    -------
    dictionary
        Dictionary object of the given case which is to be used in
        other API functions
    """

    case_dict = read_case_file(casefile)

    return case_dict


def read_case_from_mat_file(casefile):
    """
    Helper function for users
    Used to construct a case dictionary which can be used as input
    for other API functions

    Parameters
    ----------
    casefile :
        Name of and possibly full path to case file given as .mat file
        The .mat file should be in standard MATPOWER notation

    Returns
    -------
    dictionary
        Dictionary object of the given case which is to be used in
        other API functions
    """

    case_dict = read_case_file_mat(casefile)

    return case_dict


def turn_solution_into_mat_file(solution, matfilename=""):
    """
    Writes a .mat file out of an OPF solution dictionary

    Parameters
    ----------
    solution : dictionary
        OPF solution dictionary
    matfilename : string, optional
        Name of .mat file where to write the solution data
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
        "Generating .mat file %s out of given OPF solution dictionary." % matfilename
    )

    # Generate .mat file out of solution
    turn_opf_dict_into_mat_file(solution, matfilename)

    # Remove and close all logging handlers
    remove_and_close_handlers(logger, handlers)


def read_coords_from_csv_file(coordsfile):
    """
    Helper function for users
    Used to construct a coordinate dictionary which can be used as input
    for other API functions

    Parameters
    ----------
    coordsfile :
        Name of and possibly full path to case file given as .csv file

    Returns
    -------
    dictionary
        Dictionary object of the given coordinates which is to be used in
        other API functions
    """

    coord_dict = read_coords_file_csv(coordsfile)

    return coord_dict
