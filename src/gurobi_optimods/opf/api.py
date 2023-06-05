import logging

import gurobipy as gp

from gurobi_optimods.opf.grbcasereader import (
    read_case,
    read_case_file_mat,
    turn_opf_dict_into_mat_file,
)

from gurobi_optimods.opf.grbfile import (
    initialize_data_dict,
    construct_settings_dict,
    read_optimization_settings,
    read_file_csv,
    grbmap_volts_from_dict,
)
from gurobi_optimods.opf.grbformulator import (
    construct_and_solve_model,
    compute_violations_from_voltages,
)
from gurobi_optimods.opf.utils import initialize_logger, remove_and_close_handlers


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
    minactivebranches=0.9,
    useactivelossineqs=False,
    additional_settings=dict(),
):
    """
    Constructs an OPF model from given data and solves it with Gurobi.
    Returns a result dictionary following MATPOWER notation.
    The additional and possibly altered fields are

    - ``result["success"]`` 1 if at least one feasible solution has been found, 0 otherwise
    - ``result["et"]`` time spent for optimization
    - ``result["f"]`` solution objective value
    - ``result["bus"][i]["Vm"]`` for voltage magnitude value at bus `i`
    - ``result["bus"][i]["Va"]`` for voltage angle value at bus `i`
    - ``result["bus"][i]["mu"]`` for shadow prices of balance constraints at bus `i` (only available for DC without branchswitching)
    - ``result["gen"][i]["Pg"]`` for real power injection at generator `i`
    - ``result["gen"][i]["Qg"]`` for reactive power injection at generator `i`
    - ``result["branch"][i]["Pf"]`` for real power injected into "from" end of branch at branch `i`
    - ``result["branch"][i]["Pt"]`` for real power injected into "to" end of branch at branch `i`
    - ``result["branch"][i]["Qf"]`` for reactive power injected into "from" end of branch at branch `i` (AC only)
    - ``result["branch"][i]["Qt"]`` for reactive power injected into "from" end of branch at branch `i` (AC only)
    - ``result["branch"][i]["switching"]`` states whether a branch `i` is turned on or off in the final solution

    :param case: Dictionary holding case data
    :type case: dict
    :param logfile: Name of log file, defaults to ""
    :type logfile: str, optional
    :param opftype: String telling the desired OPF model type. Available are `AC`, `DC`, `IV`, defaults to `AC`
    :type opftype: str, optional
    :param polar: Controls whether polar formulation should be used, defaults to `False`. Only affects `AC` formulation
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

                            - 0 = don't use discrete variables (all branches are on)
                            - 1 = use binary variables to constrain bounds of (re)active power injections
                            - 2 = use binary variables and multiply them with the (re)active power injections (only available for AC)

                            Usually, setting 1 works better than 2. Defaults to 0
    :type branchswitching: int, optional
    :param usemipstart: Controls whether a pre-defined MIPStart should be used. Has only an effect if
                        branchswitching > 0. Defaults to `True`
    :type usemipstart: bool, optional
    :param minactivebranches: Controls the minimum number of branches that has to be turned on when branchswitching is active, i.e.,
                              the minimum number of turned on branches is equal to ``numbranches * minactivebranches``. Defaults to
                              0.95, i.e., at least 95% of branches have to be turned on
    :type minactivebranches: float, optional
    :param useactivelossineqs: Controls whether active loss constraints are used. These are linear outer approximation of the JABR
                              constraints. Usually, they provide a very good lower bound while still being linear.
                              Defaults to `False`.
    :type useactivelossineqs: bool, optional
    :param additional_settings: Dictionary holding additional settings. Additional settings are:

                                - ``lpfilename`` which if evaluated to a non-empty string, makes Gurobi write a `.lp`
                                  file holding the generated model
                                - ``gurobiparamfile`` which if evaluated to a non-empty string, makes Gurobi read in
                                  a parameter file.

                                Defaults to an empty dictionary
    :type additional_settings: dict, optional

    :return: Case dictionary following MATPOWER notation with additional result fields
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
        minactivebranches,
        useactivelossineqs,
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


def compute_violations_from_given_voltages(case, voltages, polar=False):
    """
    Constructs an OPF model from given case data and computes a violation dictionary
    out of given voltage values.

    If a voltage solution is present, e.g., from an external computation,
    this function can be used to check whether the given voltage solution is indeed
    feasible for the AC optimization model. Usually, if the violations are not too
    big, one can use the voltage solution for further calculations.

    Returns a violation dictionary following MATPOWER notation which holds
    all case data and additional violation fields. The additional fields are

    - ``violation["bus"][i]["Vmviol"]`` for voltage magnitude violation at bus `i`
    - ``violation["bus"][i]["Pviol"]`` for real injection violation at bus `i`
    - ``violation["bus"][i]["Qviol"]`` for reactive injection violation at bus `i`
    - ``violation["branch"][i]["limitviol"]`` for limit violation at branch `i`

    The violation dictionary can be used to generate a violations figure with the
    ``generate_opf_violations_figure`` function

    :param case: Dictionary holding case data
    :type case: dict
    :param voltages: Dictionary holding bus input voltage data
    :type voltages: dict
    :param polar: Controls whether polar formulation should be used when checking violations, defaults to `False`
    :type polar: bool, optional

    :return: Case dictionary following MATPOWER notation with additional violations fields
    :rtype: dict
    """

    # Initialize output and file handler and start logging
    logger, handlers = initialize_logger("OpfLogger")

    # Initialize fixed settings dictionary
    # We need the settings to construct a correct model
    settings = construct_settings_dict(
        opftype="AC",
        polar=polar,
        useef=True,
        usejabr=False,
        ivtype="aggressive",
        branchswitching=0,
        usemipstart=False,
        useactivelossineqs=False,
        minactivebranches=0.95,
        additional_settings={},
    )

    # Initilize data dictionary
    alldata = initialize_data_dict("")

    # Read settings file/dict and save them into the alldata dict
    read_optimization_settings(alldata, settings)

    # Read case file/dict and populate the alldata dictionary
    read_case(alldata, case)

    # Map given voltage data to network data
    grbmap_volts_from_dict(alldata, voltages)

    # Compute model violations based on user input voltages
    violations = compute_violations_from_voltages(alldata)

    # Remove and close all logging handlers
    remove_and_close_handlers(logger, handlers)

    return violations


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


def turn_result_into_mat_file(result, matfilename=""):
    """
    Writes a `.mat` data file in MATPOWER notation out of an OPF result dictionary.
    It accepts dictionaries in the format returned by the ``solve_opf_model`` and
    ``compute_violations_from_given_voltages`` functions

    :param result: OPF result dictionary
    :type result: dict
    :param matfilename: Name of `.mat` file where to write the result data,
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
        f"Generating .mat file {matfilename} out of given OPF result dictionary."
    )

    # Generate .mat file out of result
    turn_opf_dict_into_mat_file(result, matfilename)

    # Remove and close all logging handlers
    remove_and_close_handlers(logger, handlers)


def read_coords_from_csv_file(coordsfile):
    """
    Helper function for users.
    Constructs a coordinate dictionary which can be used as input
    for the ``generate_opf_solution_figure`` function.

    :param coordsfile: Name of and possibly full path to bus coordinates file given as `.csv` file.
                       The `.csv` file has to consist of the following columns in the given order:
                       ``index(starting with 0), busID, busname, latitude, longitude``
    :type coordsfile: str

    :return: Dictionary of the given coordinates which can be used in
             the ``generate_opf_solution_figure`` function
    :rtype: dict

    .. note::
        The first row of the csv file is ignored.
    """

    coord_dict = read_file_csv(coordsfile, "coordinates")

    return coord_dict


def read_voltages_from_csv_file(voltsfile):
    """
    Helper function for users.
    Constructs a bus input voltage dictionary which can be used as input
    for the ``check_voltage_solution_violations`` function

    :param voltsfile: Name of and possibly full path to voltage input file given as `.csv` file.
                       The `.csv` file has to consist of the following columns in the given order:
                       ``index(starting with 0), busID, busname, latitude, longitude``
    :type voltsfile: str

    :return: Dictionary of the given coordinates which can be used in
             the ``check_voltage_solution_violations`` function
    :rtype: dict

    .. note::
        The first row of the csv file is ignored.
    """

    volts_dict = read_file_csv(voltsfile, "voltages")

    return volts_dict
