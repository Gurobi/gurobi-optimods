"""
Contains the actual mods / public API.

- solve_opf: OPF solver
- compute_violations: voltage solution checker
"""

import logging

from gurobi_optimods.opf import converters, grbformulator, violations
from gurobi_optimods.utils import optimod

logger = logging.getLogger(__name__)


@optimod()
def solve_opf(
    case,
    opftype="AC",
    branch_switching=False,
    min_active_branches=0.9,
    use_mip_start=False,
    *,
    create_env,
):
    """
    Constructs an OPF model from given data and solves it with Gurobi. Returns a
    result dictionary following MATPOWER notation. The additional and possibly
    altered fields are

    - ``result["success"]`` 1 if at least one feasible solution has been found,
      0 otherwise
    - ``result["et"]`` time spent for optimization
    - ``result["f"]`` solution objective value
    - ``result["bus"][i]["Vm"]`` for voltage magnitude value at bus `i`
    - ``result["bus"][i]["Va"]`` for voltage angle value at bus `i`
    - ``result["bus"][i]["mu"]`` for shadow prices of balance constraints at bus
      `i` (only available for DC without branchswitching)
    - ``result["gen"][i]["Pg"]`` for real power injection at generator `i`
    - ``result["gen"][i]["Qg"]`` for reactive power injection at generator `i`
    - ``result["branch"][i]["Pf"]`` for real power injected into "from" end of
      branch at branch `i`
    - ``result["branch"][i]["Pt"]`` for real power injected into "to" end of
      branch at branch `i`
    - ``result["branch"][i]["Qf"]`` for reactive power injected into "from" end
      of branch at branch `i` (AC only)
    - ``result["branch"][i]["Qt"]`` for reactive power injected into "from" end
      of branch at branch `i` (AC only)
    - ``result["branch"][i]["switching"]`` states whether a branch `i` is turned
      on or off in the final solution

    Parameters
    ----------
    case : dict
        Dictionary holding case data
    opftype : str
        Desired OPF model type. One of ``AC``, ``AC_Relax``, or ``DC``
    branch_switching : bool, optional
        If set to True, enables branch switching.
    min_active_branches : float, optional
        Defines the minimum number of branches that must be turned on when
        branch switching is active, i.e. the minimum number of turned on
        branches is equal to ``numbranches * min_active_branches``. Has no
        effect if ``branchswitching`` is set to False.
    use_mip_start : bool, optional
        (Advanced) If set to True, try various MIP starts for branch switching
        models. Has no effect if ``branchswitching`` is set to False. For DC
        models, this setting is ignored, and the mip start is always used.

    Returns
    -------
    dict
        Case dictionary following MATPOWER notation with additional result
        fields
    """

    # Exact cartesian AC (force jabr for performance reasons)
    if opftype.lower() == "ac":
        opftype = "ac"
        useef = True
        usejabr = True
        default_solver_params = {"MIPGap": 1e-3, "OptimalityTol": 1e-3}
    # AC relaxation using the JABR inequality
    elif opftype.lower() == "acrelax":
        opftype = "ac"
        useef = False
        usejabr = True
        default_solver_params = {"MIPGap": 1e-3, "OptimalityTol": 1e-3}
    # DC linear approximation (ef & jabr are irrelevant)
    elif opftype.lower() == "dc":
        opftype = "dc"
        useef = False
        usejabr = False
        default_solver_params = {"MIPGap": 1e-4, "OptimalityTol": 1e-4}
    else:
        raise ValueError(f"Unknown opftype '{opftype}'")

    with create_env(params=default_solver_params) as env:
        return _solve_opf_model_internal(
            env,
            case,
            opftype=opftype,
            useef=useef,
            usejabr=usejabr,
            branchswitching=branch_switching,
            usemipstart=use_mip_start,
            minactivebranches=min_active_branches,
            polar=False,
            ivtype="aggressive",
            useactivelossineqs=False,
        )


def _solve_opf_model_internal(
    env,
    case,
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
    # Initialize settings dictionary
    settings = converters.build_internal_settings(
        opftype,
        polar,
        useef,
        usejabr,
        ivtype,
        branchswitching,
        usemipstart,
        minactivebranches,
        useactivelossineqs,
    )

    # Populate the alldata dictionary with case data and settings
    alldata = converters.convert_case_to_internal_format(case)
    alldata.update(settings)

    # Construct and solve model using given case data and user settings
    solution = grbformulator.construct_and_solve_model(env, alldata)

    # TODO solution data is populated into 'alldata' then extracted by another
    # function. It may make more sense to have construct_and_solve_model return
    # nothing, then extract the solution as a separate process. This is more
    # consistent with the way the rest of the code works (essentially dumping
    # everything into 'alldata').

    return solution


@optimod()
def compute_violations(case, voltages, polar=False, *, create_env):
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

    Parameters
    ----------
    case : dict
        Dictionary holding case data
    voltages : dict
        Dictionary holding bus input voltage data
    polar: bool, optional
        If True, use polar formulation when checking violations, defaults to False

    Returns
    -------
    dict
        Case dictionary following MATPOWER notation with additional violations fields

    """

    # Initialize fixed settings dictionary
    # We need the settings to construct a correct model
    settings = converters.build_internal_settings(
        opftype="AC",
        polar=polar,
        useef=True,
        usejabr=False,
        ivtype="aggressive",
        branchswitching=0,
        usemipstart=False,
        useactivelossineqs=False,
        minactivebranches=0.95,
    )

    # Populate the alldata dictionary with case data and settings
    alldata = converters.convert_case_to_internal_format(case)
    alldata.update(settings)

    # Map given voltage data to network data
    converters.grbmap_volts_from_dict(alldata, voltages)

    # Compute model violations based on user input voltages
    with create_env() as env:
        return violations.compute_violations_from_voltages(env, alldata)
