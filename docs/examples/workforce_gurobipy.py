import pandas as pd

import gurobipy as gp


availability = pd.read_csv("data/availability.csv").assign(
    Shift=lambda df: pd.to_datetime(df["Shift"])
)
shift_requirements = (
    pd.read_csv("data/shift_requirements.csv")
    .assign(Shift=lambda df: pd.to_datetime(df["Shift"]))
    .set_index("Shift")["Required"]
)
pay_rates = pd.read_csv("data/pay_rates.csv").set_index("Worker")["PayRate"]

m = gp.Model()

x = m.addMVar(availability.shape[0], ub=1)

# Objective:
for worker, worker_shifts in availability.groupby("Worker"):
    x[worker_shifts.index].Obj = pay_rates.loc[worker]

# Constraint: enough workers per shifts
for shift, shift_workers in availability.groupby("Shift"):
    m.addConstr(x[shift_workers.index].sum() == shift_requirements.loc[shift])

m.optimize()

# TODO fix data so shift ordering is sensible; use dates?

# Use solution to filter selected shifts.
assigned_shifts = availability[pd.Series(index=availability.index, data=x.X > 0.9)]
