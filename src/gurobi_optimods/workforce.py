from typing import Optional

import gurobipy as gp
import gurobipy_pandas as gppd
import pandas as pd
from gurobipy import GRB


def solve_workforce_scheduling(
    availability: pd.DataFrame,
    shift_requirements: pd.DataFrame,
    worker_data: Optional[pd.DataFrame] = None,
    rolling_window=None,
    rolling_limit=None,
) -> pd.DataFrame:
    """Solve a workforce scheduling model.

    :param availability: Dataframe with columns 'Worker' and 'Shift' defining
        all allowable worker-shift combinations
    :type availability: :class:`pd.DataFrame`
    :param shift_requirements: Dataframe with columns 'Shift' and 'Required'
        specifying the number of staff required for every shift
    :type shift_requirements: :class:`pd.DataFrame`
    """
    with gp.Env() as env, gp.Model(env=env) as m:
        m.ModelSense = GRB.MAXIMIZE
        assignments = availability.set_index(["Worker", "Shift"]).gppd.add_vars(
            m, obj="Preference", vtype=GRB.BINARY, name="assign"
        )
        gppd.add_constrs(
            m,
            assignments.groupby("Shift")["assign"].sum(),
            GRB.EQUAL,
            shift_requirements.set_index("Shift")["Required"],
            name="requirements",
        )

        if worker_data is not None:
            gppd.add_constrs(
                m,
                assignments.groupby("Worker")["assign"].sum(),
                GRB.LESS_EQUAL,
                worker_data.set_index("Worker")["MaxShifts"],
                name="max_shifts",
            )
            gppd.add_constrs(
                m,
                assignments.groupby("Worker")["assign"].sum(),
                GRB.GREATER_EQUAL,
                worker_data.set_index("Worker")["MinShifts"],
                name="min_shifts",
            )

        if rolling_window:
            # This is uglier than it should be, but pandas rolling() only
            # handles numeric data
            for _, df in assignments.reset_index().groupby("Worker"):
                df = df.set_index("Shift")["assign"]
                for entry in df.index:
                    # Hack! Surely there is a open/closed interval available?
                    expr = df.loc[
                        entry : entry + rolling_window - pd.Timedelta(seconds=1)
                    ].sum()
                    m.addConstr(expr <= rolling_limit)

            # TODO need another option to specify a rolling roster (wraparound)

        m.optimize()
        return (
            assignments.assign(assign=lambda df: df["assign"].gppd.X)
            .query("assign > 0.9")
            .drop(columns=["assign"])
            .reset_index()
        )
