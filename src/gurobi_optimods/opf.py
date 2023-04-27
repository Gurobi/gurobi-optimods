import logging

import gurobipy as gp

from .src_opf.grbcasereader import (
    read_case,
    read_case_file_mat,
    turn_opf_dict_into_mat_file,
)

from .src_opf.grbfile import (
    initialize_data_dict,
    construct_settings_dict,
    read_optimization_settings,
    read_graphics_settings,
    read_coords_file_csv,
    grbmap_coords_from_dict,
    grbread_graphattrs,
)
from .src_opf.grbformulator import construct_and_solve_model
from .src_opf.grbgraphical import generate_solution_figure
from .src_opf.utils import initialize_logger, remove_and_close_handlers


def solve_opf_model(
    case,
    logfile="",
    opftype="AC",
    polar=False,
    useef=True,
    usejabr=True,
    ivtype="aggressive",
    branchswitching=0,
    usemipstart=True,
    additional_settings=dict(),
):
    """
    Constructs an OPF model from given data and solves it with Gurobi.
    Returns a result dictionary following MATPOWER notation

    :param case: Dictionary holding case data
    :type case: dict
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

    :return: A result dictionary following MATPOWER notation, e.g., the objective value
             is at ``solution["f"]`` if ``solution["success"] = 1``
    :rtype: dict
    """

    # Initialize output and file handler and start logging
    logger, handlers = initialize_logger("OpfLogger", logfile, True)

    # Initialize settings dictionary
    settings = construct_settings_dict(
        opftype,
        polar,
        useef,
        usejabr,
        ivtype,
        branchswitching,
        usemipstart,
        additional_settings,
    )

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


def generate_opf_solution_figure(
    case, coords, solution, graphattrsfile="", voltsfile=""
):
    """
    Reads the given case and returns a plotly figure object.
    Ideally the solution has been computed by the `solve_opf_model` function

    :param case: Dictionary holding case data
    :type case: dict
    :param coords: Dictionary holding bus coordinates
    :type coords: dict
    :param solution: Dictionary holding solution data following the MATPOWER notation as returned
                      by the solve_opf_model function
    :type solution: dict
    :param graphattrsfile: Name of and possibly full path to a file holding additional graph attributes,
                           defaults to "". #TODO will be very likely removed
    :type graphattrsfile: str, optional
    :param voltsfile: Name of and possibly full path to a file holding voltage solution information,
                      defaults to "". #TODO will be very likely removed
    :type voltsfile: str, optional

    :return: A plotly figure objects which can be displayed via the show() function,
             see https://plotly.com/python-api-reference/generated/plotly.graph_objects.Figure.html
    :rtype: :class:`plotly.graph_objects.Figure`
    """

    # Initialize output and file handler and start logging
    logger, handlers = initialize_logger("OpfLogger")

    # Initilize data dictionary
    alldata = initialize_data_dict()

    settings = {"voltsfilename": voltsfile, "graphattrsfilename": graphattrsfile}

    # Read settings file/dict and save them into the alldata dict
    read_graphics_settings(alldata, settings)

    # Read case file/dict and populate the alldata dictionary
    read_case(alldata, case)

    # Special settings for graphics
    alldata["graphical"] = {}
    alldata["graphical"]["numfeatures"] = 0
    if alldata["graphattrsfilename"] != "" and alldata["graphattrsfilename"] != None:
        grbread_graphattrs(alldata, alldata["graphattrsfilename"])

    # Map given coordinate data to network data
    grbmap_coords_from_dict(alldata, coords)

    # Generate a plotly figure object representing the given solution for the network
    fig = generate_solution_figure(alldata, solution)

    # Remove and close all logging handlers
    remove_and_close_handlers(logger, handlers)

    return fig


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
    :rtype: dict
    """

    case_dict = read_case_file_mat(casefile)

    return case_dict


def turn_solution_into_mat_file(solution, matfilename=""):
    """
    Writes a `.mat` data file in MATPOWER notation out of an OPF solution dictionary

    :param solution: OPF solution dictionary
    :type solution: dict
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
    Constructs a coordinate dictionary which can be used as input
    for the `generate_opf_solution_figure` function

    :param coordsfile: Name of and possibly full path to case file given as `.csv` file
    :type coordsfile: str

    :return: Dictionary of the given coordinates which can be used in
             the ``generate_opf_solution_figure`` function
    :rtype: dict
    """

    coord_dict = read_coords_file_csv(coordsfile)

    return coord_dict
