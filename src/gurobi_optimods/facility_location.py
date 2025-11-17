"""
Facility Location
-----------------
"""

import logging

import gurobipy as gp
import gurobipy_pandas as gppd
from gurobipy import GRB

from gurobi_optimods.utils import optimod

logger = logging.getLogger(__name__)


@optimod()
def solve_facility_location(
    customer_data,
    facility_data,
    transportation_cost,
    *,
    fixed_facility_count=None,
    max_facilities_per_region=None,
    create_env,
):
    """Solve a capacitated facility location problem with multiple allocations.

    Given a set of customers with known demand, a set of potential facility
    locations with fixed opening costs and capacities, and transportation costs
    between customers and facilities, determine which facilities to open and how
    to allocate customer demand to facilities to minimize total cost.

    Parameters
    ----------
    customer_data : DataFrame
        Customer information with columns:

        - ``customer``: Customer identifier
        - ``demand``: Demand quantity for this customer (must be non-negative)

    facility_data : DataFrame
        Facility information with columns:

        - ``facility``: Facility identifier
        - ``capacity``: Maximum capacity of this facility (must be non-negative)
        - ``fixed_cost``: Fixed cost to open this facility (must be non-negative)
        - ``region`` (optional): Region or category identifier for fairness constraints

    transportation_cost : DataFrame
        Transportation cost data with columns:

        - ``customer``: Customer identifier
        - ``facility``: Facility identifier
        - ``cost``: Cost per unit to transport from facility to customer (must be non-negative)

    fixed_facility_count : int, optional
        If specified, exactly this many facilities must be opened. If not specified,
        the number of facilities to open is determined by the optimization to minimize
        total cost. Must be a positive integer not exceeding the total number of facilities.

    max_facilities_per_region : int, optional
        If specified, limits the maximum number of facilities that can be opened in
        any single region. Requires the ``region`` column in ``facility_data``.
        Use case: Prevent concentration of undesirable facilities (e.g., landfills)
        in specific regions or demographics to ensure equitable distribution.

    Returns
    -------
    dict
        A dictionary containing:

        - ``solution_value``: Total cost of the optimal solution (fixed + transportation costs)
        - ``facilities_opened``: A pandas Series indexed by facility showing which facilities are opened (1.0 if open, 0.0 if closed)
        - ``assignments``: A pandas DataFrame with columns ``customer``, ``facility``, and ``assignment`` showing the quantity allocated from each facility to each customer

    Raises
    ------
    ValueError
        If input data is missing required columns, contains negative values where
        not allowed, is empty, if the problem is infeasible, if fixed_facility_count
        is invalid, or if max_facilities_per_region is improperly specified.
    """

    # Validate inputs
    _validate_inputs(
        customer_data,
        facility_data,
        transportation_cost,
        fixed_facility_count,
        max_facilities_per_region,
    )

    # Build and solve the model
    with create_env() as env, gp.Model(env=env) as model:
        # Decision variables
        # Binary variable: is facility j open?
        facility_df = facility_data.set_index("facility")
        facility_open = gppd.add_vars(
            model, facility_df, name="facility_open", vtype=GRB.BINARY
        )

        # Continuous variable: amount shipped from facility j to customer i
        transport_df = transportation_cost.set_index(["customer", "facility"])
        transport_vars = gppd.add_vars(model, transport_df, name="transport", lb=0.0)

        # Objective: minimize total cost (fixed costs + transportation costs)
        fixed_costs = (facility_df["fixed_cost"] * facility_open).sum()
        transport_costs = (transport_df["cost"] * transport_vars).sum()
        model.setObjective(fixed_costs + transport_costs, GRB.MINIMIZE)

        # Constraint: satisfy all customer demand
        customer_demand = customer_data.set_index("customer")["demand"]
        transport_by_customer = transport_vars.groupby(level="customer").sum()

        gppd.add_constrs(
            model,
            transport_by_customer,
            GRB.EQUAL,
            customer_demand,
            name="satisfy_demand",
        )

        # Constraint: respect facility capacity
        facility_capacity = facility_df["capacity"]
        transport_by_facility = transport_vars.groupby(level="facility").sum()

        gppd.add_constrs(
            model,
            transport_by_facility,
            GRB.LESS_EQUAL,
            facility_capacity,
            name="capacity_limit",
        )

        # Constraint: can only ship from open facilities (big-M constraint)
        # For each (customer, facility) pair: transport[i,j] <= demand[i] * facility_open[j]
        # This ensures that if facility j is closed, no shipment can go from j to any customer
        for idx in transport_df.index:
            customer, facility = idx
            customer_demand_value = customer_demand.loc[customer]

            model.addConstr(
                transport_vars.loc[idx]
                <= customer_demand_value * facility_open.loc[facility],
                name=f"open_facility_{customer}_{facility}",
            )

        # Constraint: fixed number of facilities (if specified)
        if fixed_facility_count is not None:
            model.addConstr(
                facility_open.sum() == fixed_facility_count,
                name="fixed_facility_count",
            )

        # Constraint: fairness limits per region (if specified)
        if max_facilities_per_region is not None:
            facility_regions = facility_data.set_index("facility")["region"]

            for region in facility_regions.unique():
                facilities_in_region = facility_regions[
                    facility_regions == region
                ].index
                model.addConstr(
                    facility_open.loc[facilities_in_region].sum()
                    <= max_facilities_per_region,
                    name=f"fairness_region_{region}",
                )

        # Solve the model
        model.optimize()

        # Check solution status
        if model.status == GRB.INFEASIBLE:
            raise ValueError(
                "The facility location problem is infeasible. "
                "This may occur if total facility capacity is insufficient to meet total demand."
            )
        elif model.status != GRB.OPTIMAL:
            raise ValueError(
                f"Optimization failed with status {model.status}. "
                "Please check your input data."
            )

        # Extract solution
        solution_value = model.ObjVal

        # Get opened facilities (using .gppd.X accessor)
        facilities_opened = facility_open.gppd.X

        # Get assignments (filter out zero assignments for cleaner output)
        assignments = transport_vars.gppd.X
        assignments = assignments[assignments > 1e-6].reset_index()
        assignments.columns = ["customer", "facility", "assignment"]

        logger.info(f"Optimal solution found with total cost: {solution_value:.2f}")
        logger.info(
            f"Facilities opened: {facilities_opened[facilities_opened > 0.5].index.tolist()}"
        )

        return {
            "solution_value": solution_value,
            "facilities_opened": facilities_opened,
            "assignments": assignments,
        }


def _validate_inputs(
    customer_data,
    facility_data,
    transportation_cost,
    fixed_facility_count=None,
    max_facilities_per_region=None,
):
    """Validate input dataframes for facility location problem.

    Parameters
    ----------
    customer_data : DataFrame
        Customer information
    facility_data : DataFrame
        Facility information
    transportation_cost : DataFrame
        Transportation cost data
    fixed_facility_count : int, optional
        Fixed number of facilities to open
    max_facilities_per_region : int, optional
        Maximum facilities allowed per region

    Raises
    ------
    ValueError
        If any validation check fails
    """

    # Check for empty dataframes
    if customer_data.empty:
        raise ValueError("customer_data cannot be empty")
    if facility_data.empty:
        raise ValueError("facility_data cannot be empty")
    if transportation_cost.empty:
        raise ValueError("transportation_cost cannot be empty")

    # Check required columns in customer_data
    required_customer_cols = ["customer", "demand"]
    missing_cols = set(required_customer_cols) - set(customer_data.columns)
    if missing_cols:
        raise ValueError(f"customer_data is missing required columns: {missing_cols}")

    # Check required columns in facility_data
    required_facility_cols = ["facility", "capacity", "fixed_cost"]
    missing_cols = set(required_facility_cols) - set(facility_data.columns)
    if missing_cols:
        raise ValueError(f"facility_data is missing required columns: {missing_cols}")

    # Check required columns in transportation_cost
    required_transport_cols = ["customer", "facility", "cost"]
    missing_cols = set(required_transport_cols) - set(transportation_cost.columns)
    if missing_cols:
        raise ValueError(
            f"transportation_cost is missing required columns: {missing_cols}"
        )

    # Check for negative demands
    if (customer_data["demand"] < 0).any():
        raise ValueError("All demand values must be non-negative")

    # Check for negative capacities
    if (facility_data["capacity"] < 0).any():
        raise ValueError("All capacity values must be non-negative")

    # Check for negative fixed costs
    if (facility_data["fixed_cost"] < 0).any():
        raise ValueError("All fixed_cost values must be non-negative")

    # Check for negative transportation costs
    if (transportation_cost["cost"] < 0).any():
        raise ValueError("All transportation cost values must be non-negative")

    # Validate fixed_facility_count if provided
    if fixed_facility_count is not None:
        if not isinstance(fixed_facility_count, int):
            raise ValueError("fixed_facility_count must be an integer")
        if fixed_facility_count <= 0:
            raise ValueError("fixed_facility_count must be a positive integer")
        num_facilities = len(facility_data)
        if fixed_facility_count > num_facilities:
            raise ValueError(
                f"fixed_facility_count ({fixed_facility_count}) exceeds "
                f"the number of available facilities ({num_facilities})"
            )

    # Validate max_facilities_per_region if provided
    if max_facilities_per_region is not None:
        if not isinstance(max_facilities_per_region, int):
            raise ValueError("max_facilities_per_region must be an integer")
        if max_facilities_per_region <= 0:
            raise ValueError("max_facilities_per_region must be a positive integer")

        # Check that facility_data has region column
        if "region" not in facility_data.columns:
            raise ValueError(
                "facility_data must contain 'region' column when "
                "max_facilities_per_region is specified"
            )
