import pandas as pd

import gurobipy as gp


availability = pd.read_feather("examples/data/availability.feather")
shift_requirements = (
    pd.read_feather("examples/data/shift_requirements.feather")
    .set_index("Shift")["Required"]
)
pay_rates = (
    pd.read_feather("examples/data/pay_rates.feather")
    .set_index("Worker")["PayRate"]
)

m = gp.Model()

x = m.addMVar(availability.shape[0], ub=1)

# Objective:
for worker, worker_shifts in availability.groupby("Worker"):
    x[worker_shifts.index].Obj = pay_rates.loc[worker]

# Constraint: enough workers per shifts
for shift, shift_workers in availability.groupby("Shift"):
    m.addConstr(x[shift_workers.index].sum() == shift_requirements.loc[shift])

m.optimize()

# Use solution to filter selected shifts.
assigned_shifts = availability[pd.Series(index=availability.index, data=x.X > 0.9)]
