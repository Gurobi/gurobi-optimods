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
    worker_limits: Optional[pd.DataFrame] = None,
    rolling_window=None,
    rolling_limit=None,
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
    """
    with create_env() as env, gp.Model(env=env) as m:
        m.ModelSense = GRB.MAXIMIZE
        assignments = preferences.set_index(["Worker", "Shift"]).gppd.add_vars(
            m, obj="Preference", vtype=GRB.BINARY, name="assign"
        )
        gppd.add_constrs(
            m,
            assignments.groupby("Shift")["assign"].sum(),
            GRB.EQUAL,
            shift_requirements.set_index("Shift")["Required"],
            name="requirements",
        )

        if worker_limits is not None:
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

        if rolling_window:
            # This is uglier than it should be, but pandas rolling() only
            # handles numeric data
            for worker, df in assignments.reset_index().groupby("Worker"):
                df = df.set_index("Shift")["assign"]
                for entry in df.index:
                    # Hack! Surely there is a open/closed interval available?
                    expr = df.loc[
                        entry : entry + rolling_window - pd.Timedelta(seconds=1)
                    ].sum()
                    m.addConstr(
                        expr <= rolling_limit, name=f"rolling[{worker},{entry}]"
                    )

            # TODO need another option to specify a rolling roster (wraparound)

        m.optimize()

        if m.Status == GRB.INFEASIBLE:
            raise ValueError("Infeasible roster")

        return (
            assignments.assign(assign=lambda df: df["assign"].gppd.X)
            .query("assign > 0.9")
            .drop(columns=["assign"])
            .reset_index()
        )
