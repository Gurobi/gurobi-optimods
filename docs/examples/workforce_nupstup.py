import pandas as pd

from nupstup.workforce import solve_workforce_scheduling


availability = pd.read_csv("data/availability.csv").assign(
    Shift=lambda df: pd.to_datetime(df["Shift"])
)
shift_requirements = (
    pd.read_csv("data/shift_requirements.csv")
    .assign(Shift=lambda df: pd.to_datetime(df["Shift"]))
    .set_index("Shift")["Required"]
)
pay_rates = pd.read_csv("data/pay_rates.csv").set_index("Worker")["PayRate"]

# Get results directly from standard data input format.
# TODO: should be a Model or Solver class to provide options
# (e.g. LogFile, LogToConsole)?
assigned_shifts = solve_workforce_scheduling(
    availability=availability,
    shift_requirements=shift_requirements,
    pay_rates=pay_rates,
)
