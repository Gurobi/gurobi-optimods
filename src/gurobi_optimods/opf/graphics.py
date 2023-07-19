"""
Contains the plotting API for OPF.

- plot_solution: plot case solution using plotly, highlighting switched off branches
- plot_violations: plot case solution, highlighting violations

TODO: aren't these two basically the same? 'violations' is just a solution with
some extra fields. plot_solution(show_violations=True) might be better?

TODO: python users expect functions called plot_ e.g. plot_solution,
plot_violations
"""


from gurobi_optimods.opf import converters, grbgraphical


def plot_solution(case, coords, solution):
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

    # Populate the alldata dictionary with case data
    alldata = converters.convert_case_to_internal_format(case)

    # Special settings for graphics
    alldata["graphical"] = {}
    alldata["graphical"]["numfeatures"] = 0

    # Map given coordinate data to network data
    converters.grbmap_coords_from_dict(alldata, coords)

    # Generate a plotly figure object representing the given solution for the network
    fig = grbgraphical.generate_solution_figure(alldata, solution)

    return fig


def plot_violations(case, coords, violations):
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

    # Populate the alldata dictionary with case data
    alldata = converters.convert_case_to_internal_format(case)

    # Special settings for graphics
    alldata["graphical"] = {}
    alldata["graphical"]["numfeatures"] = 0

    # Map given coordinate data to network data
    converters.grbmap_coords_from_dict(alldata, coords)

    # Generate a plotly figure object representing the given violations for the network
    fig = grbgraphical.generate_violations_figure(alldata, violations)

    return fig
