import pandas as pd

from nupstup.workforce import solve_workforce_scheduling


availability = pd.read_csv("data/availability.csv")
shift_requirements = pd.read_csv("data/shiftReq.csv", index_col=[0])
hourly_rates = pd.read_csv("data/workerpay.csv", index_col=[0])

# Get results directly from standard data input format.
assigned_shifts = solve_workforce_scheduling(
    availability=availability,
    shift_requirements=shift_requirements,
    hourly_rates=hourly_rates,
)
