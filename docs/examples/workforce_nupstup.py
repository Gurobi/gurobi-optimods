import pandas as pd

from nupstup.workforce import solve_workforce_scheduling


# Load example data.
availability = pd.read_feather("examples/data/availability.feather")
shift_requirements = pd.read_feather("examples/data/shift_requirements.feather")
pay_rates = pd.read_feather("examples/data/pay_rates.feather")

# Get winning results.
assigned_shifts = solve_workforce_scheduling(
    availability=availability,
    shift_requirements=shift_requirements,
    pay_rates=pay_rates,
)
