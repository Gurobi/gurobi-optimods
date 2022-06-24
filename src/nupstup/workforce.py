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
    A = sp.coo_matrix((data, (row, col)))
    b = shift_req.loc[categorical.cat.categories].values
    m.addMConstr(A, x, gp.GRB.EQUAL, b)

    return x


def solve_workforce_scheduling(
    availability: pd.DataFrame, shift_requirements: pd.Series, pay_rates: pd.Series
) -> pd.DataFrame:
    """Solve a workforce scheduling model"""
    m = gp.Model()
    x = workforce_mconstr(m, availability, shift_requirements, pay_rates)
    m.optimize()
    return availability[pd.Series(index=availability.index, data=x.X > 0.9)].reset_index(drop=True)
