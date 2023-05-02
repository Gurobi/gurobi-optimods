from gurobi_optimods.opf.grbgraphical import (
    generate_solution_figure,
    generate_violations_figure,
)
from gurobi_optimods.opf.utils import initialize_logger, remove_and_close_handlers
from gurobi_optimods.opf.grbfile import (
    initialize_data_dict,
    grbmap_coords_from_dict,
)
from gurobi_optimods.opf.grbcasereader import read_case


def generate_opf_solution_figure(case, coords, solution):
    """
    Reads the given case and returns a plotly figure object.
    Ideally the solution has been computed by the ``solve_opf_model`` function

    :param case: Dictionary holding case data
    :type case: dict
    :param coords: Dictionary holding bus coordinates
    :type coords: dict
    :param solution: Dictionary holding solution data following the MATPOWER notation as returned
                      by the ``solve_opf_model`` function
    :type solution: dict

    :return: A plotly figure objects which can be displayed via the show() function,
             see https://plotly.com/python-api-reference/generated/plotly.graph_objects.Figure.html
    :rtype: :class:`plotly.graph_objects.Figure`
    """

    # Initialize output and file handler and start logging
    logger, handlers = initialize_logger("OpfLogger")

    # Initilize data dictionary
    alldata = initialize_data_dict()

    # Read case file/dict and populate the alldata dictionary
    read_case(alldata, case)

    # Special settings for graphics
    alldata["graphical"] = {}
    alldata["graphical"]["numfeatures"] = 0

    # Map given coordinate data to network data
    grbmap_coords_from_dict(alldata, coords)

    # Generate a plotly figure object representing the given solution for the network
    fig = generate_solution_figure(alldata, solution)

    # Remove and close all logging handlers
    remove_and_close_handlers(logger, handlers)

    return fig


def generate_opf_violations_figure(case, coords, violations):
    """
    Reads the given case and returns a plotly figure object of provided violations.
    Ideally the violations have been computed by the ``compute_violations_from_given_voltages`` function

    :param case: Dictionary holding case data
    :type case: dict
    :param coords: Dictionary holding bus coordinates
    :type coords: dict
    :param violations: Dictionary holding case data following the MATPOWER notation with additional
                       violations fields as returned by the ``compute_violations_from_given_voltages`` function
    :type violations: dict

    :return: A plotly figure objects which can be displayed via the show() function,
             see https://plotly.com/python-api-reference/generated/plotly.graph_objects.Figure.html
    :rtype: :class:`plotly.graph_objects.Figure`
    """

    # Initialize output and file handler and start logging
    logger, handlers = initialize_logger("OpfLogger")

    # Initilize data dictionary
    alldata = initialize_data_dict()

    # Read case file/dict and populate the alldata dictionary
    read_case(alldata, case)

    # Special settings for graphics
    alldata["graphical"] = {}
    alldata["graphical"]["numfeatures"] = 0

    # Map given coordinate data to network data
    grbmap_coords_from_dict(alldata, coords)

    # Generate a plotly figure object representing the given violations for the network
    fig = generate_violations_figure(alldata, violations)

    # Remove and close all logging handlers
    remove_and_close_handlers(logger, handlers)

    return fig
