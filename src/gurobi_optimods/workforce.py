import logging
from typing import Optional

import gurobipy as gp
import gurobipy_pandas as gppd
import pandas as pd
from gurobipy import GRB

from gurobi_optimods.utils import optimod

logger = logging.getLogger(__name__)


@optimod()
def solve_workforce_scheduling(
    preferences: pd.DataFrame,
    shift_requirements: pd.DataFrame,
    worker_limits: pd.DataFrame,
    limit_window: Optional[pd.Timedelta] = None,
    *,
    create_env,
) -> pd.DataFrame:
    """Solve a workforce scheduling model.

    :param preferences: Dataframe with columns 'Worker' and 'Shift' defining
        all allowable worker-shift combinations
    :type preferences: :class:`pd.DataFrame`
    :param shift_requirements: Dataframe with columns 'Shift' and 'Required'
        specifying the number of staff required for every shift
    :type shift_requirements: :class:`pd.DataFrame`
    :param worker_limits: Dataframe with columns 'Worker', 'MinShifts', and
        'MaxShifts' specifying the maximum and minimum number of shifts each
        worker may be assigned
    :type worker_limits: :class:`pd.DataFrame`
    :param limit_window: (Optional) Rolling window
    :type limit_window: :class:`pd.Timedelta`
    :return: Assigned shifts as a subset of the preferences dataframe
    :rtype: :class:`pd.DataFrame`
    """
    with create_env() as env, gp.Model(env=env) as m:

        # Create binary variables for all valid shift assignments and
        # create preference maximization objective
        m.ModelSense = GRB.MAXIMIZE
        assignments = preferences.set_index(["Worker", "Shift"]).gppd.add_vars(
            m, obj="Preference", vtype=GRB.BINARY, name="assign"
        )

        # Enforce shift coverage requirements
        gppd.add_constrs(
            m,
            assignments.groupby("Shift")["assign"].sum(),
            GRB.EQUAL,
            shift_requirements.set_index("Shift")["Required"],
            name="requirements",
        )

        if limit_window:
            # If limit_window is specified, min/max shift limits are interpreted
            # as limits on rolling windows of the roster, where windows length
            # is dictated by the limit_window value
            for worker, df in assignments.reset_index().groupby("Worker"):
                df = df.set_index("Shift")["assign"]
                for entry in df.index:
                    # TODO hacky! Surely there is a open/closed interval available?
                    expr = df.loc[
                        entry : entry + limit_window - pd.Timedelta(seconds=1)
                    ].sum()
                    m.addConstr(
                        expr
                        <= worker_limits.set_index("Worker").loc[worker, "MaxShifts"],
                        name=f"rolling[{worker},{entry}]",
                    )
                    # TODO test and implement lower limit

        else:
            # If limit_window is not specified, min/max shift limits are
            # interpreted as limits on the roster as a whole
            gppd.add_constrs(
                m,
                assignments.groupby("Worker")["assign"].sum(),
                GRB.LESS_EQUAL,
                worker_limits.set_index("Worker")["MaxShifts"],
                name="max_shifts",
            )
            gppd.add_constrs(
                m,
                assignments.groupby("Worker")["assign"].sum(),
                GRB.GREATER_EQUAL,
                worker_limits.set_index("Worker")["MinShifts"],
                name="min_shifts",
            )

        # Solve the model and return the shift assignments as a subset of the
        # input preferences dataframe. Raise an exception if a feasible schedule
        # does not exist.

        m.optimize()
        if m.Status == GRB.INFEASIBLE:
            raise ValueError("Infeasible roster")

        return (
            assignments.assign(assign=lambda df: df["assign"].gppd.X)
            .query("assign > 0.9")
            .drop(columns=["assign"])
            .reset_index()
        )
