"""
Contains the plotting API for OPF.

- solution_plot: plot case solution using plotly, highlighting switched off branches
- violation_plot: plot case solution, highlighting violations

"""


from gurobi_optimods.opf import converters, grbgraphical


def solution_plot(case, coords, solution):
    """
    Reads the given case and returns a plotly figure object. Ideally the
    solution has been computed by the ``solve_opf`` function

    Parameters
    ----------
    case : dict
        Dictionary holding case data
    coords : dict
        Dictionary holding bus coordinates
    solution: dict
        Dictionary holding solution data following the MATPOWER notation as
        returned by the ``solve_opf`` function

    Returns
    -------
    plotly.graph_objects.Figure
        A plotly figure object displaying the solution. The plot can be
        displaged by calling ``figure.show()``.
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


def violation_plot(case, coords, violations):
    """
    Reads the given case and returns a plotly figure object of provided
    violations. Ideally the violations have been computed by the
    ``compute_violations`` function

    Parameters
    ----------
    case : dict
        Dictionary holding case data
    coords : dict
        Dictionary holding bus coordinates
    violations : dict
        Dictionary holding case data following the MATPOWER notation with
        additional violations fields as returned by the ``compute_violations``
        function

    Returns
    -------
    plotly.graph_objects.Figure
        A plotly figure object highlighting violations in the solution. The
        plot can be displaged by calling ``figure.show()``.
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
