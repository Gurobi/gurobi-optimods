import gurobipy as gp
import gurobipy_pandas as gppd
import pandas as pd
from gurobipy import GRB


def solve_workforce_scheduling(
    availability: pd.DataFrame, shift_requirements: pd.DataFrame
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
        m.optimize()
        return (
            assignments.assign(assign=lambda df: df["assign"].gppd.X)
            .query("assign > 0.9")
            .drop(columns=["assign"])
            .reset_index()
        )
