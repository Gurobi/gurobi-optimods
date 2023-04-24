import gurobipy as gp
import numpy as np
import pandas as pd
import scipy.sparse as sp


def workforce_mconstr(m, availability, shift_req, pay):

    # Directly set the objective at construction time.
    x = m.addMVar(availability.index.size, ub=1, obj=pay.loc[availability["Worker"]])

    # Use pandas to munge data into matrix form.
    categorical = availability["Shift"].astype("category")
    row = categorical.cat.codes
    col = availability.index
    data = np.ones(availability.index.size)
    A = sp.coo_array((data, (row, col)))
    b = shift_req.loc[categorical.cat.categories].values
    m.addMConstr(A, x, gp.GRB.EQUAL, b)

    return x


def solve_workforce_scheduling(
    availability: pd.DataFrame, shift_requirements: pd.Series, pay_rates: pd.Series
) -> pd.DataFrame:
    """Solve a workforce scheduling model.

    :param availability: Dataframe with columns 'Worker' and 'Shift' defining
        all allowable worker-shift combinations
    :type availability: :class:`pd.DataFrame`
    :param shift_requirements: Dataframe with columns 'Shift' and 'Required'
        specifying the number of staff required for every shift
    :type shift_requirements: :class:`pd.DataFrame`
    :param pay_rates: Dataframe with columns 'Worker' and 'PayRate' specifying
        the per-shift pay rate of every worker
    :type pay_rates: :class:`pd.DataFrame`
    """
    with gp.Env() as env, gp.Model(env=env) as m:
        x = workforce_mconstr(
            m,
            availability,
            shift_requirements.set_index("Shift")["Required"],
            pay_rates.set_index("Worker")["PayRate"],
        )
        m.optimize()
        return availability[
            pd.Series(index=availability.index, data=x.X > 0.9)
        ].reset_index(drop=True)
