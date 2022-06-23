import pandas as pd

import gurobipy as gp


availability = pd.read_csv("data/availability.csv")
shift_req = pd.read_csv("data/shiftReq.csv", index_col=[0])
pay = pd.read_csv("data/workerpay.csv", index_col=[0])

m = gp.Model()

x = m.addMVar(availability.shape[0], ub=1)

# Objective:
for worker, worker_shifts in availability.groupby("Workers"):
    x[worker_shifts.index].Obj = pay.loc[worker, "Pay"]

# Constraint: enough workers per shifts
for shift, shift_workers in availability.groupby("Shift"):
    m.addConstr(x[shift_workers.index].sum() == shift_req.loc[shift])

m.optimize()

# TODO fix data so shift ordering is sensible; use dates?

# Use solution to filter selected shifts.
assigned_shifts = availability[pd.Series(index=availability.index, data=x.X > 0.9)]
