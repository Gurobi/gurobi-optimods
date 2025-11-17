import unittest

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal, assert_series_equal

from gurobi_optimods.facility_location import solve_facility_location


class TestFacilityLocationBasic(unittest.TestCase):
    """Test basic functionality of facility location optimod."""

    def test_simple_uncapacitated_case(self):
        """Test a simple case where solution is obvious."""
        # 2 customers, 2 facilities
        # Facility 0: fixed cost 10, serves customer 0 (demand 5) at cost 1 per unit
        # Facility 1: fixed cost 10, serves customer 1 (demand 5) at cost 1 per unit
        # Optimal: open both facilities
        # Total cost = 10 + 10 + (5*1) + (5*1) = 30

        customer_data = pd.DataFrame({"customer": [0, 1], "demand": [5.0, 5.0]})

        facility_data = pd.DataFrame(
            {"facility": [0, 1], "capacity": [10.0, 10.0], "fixed_cost": [10.0, 10.0]}
        )

        transportation_cost = pd.DataFrame(
            {
                "customer": [0, 0, 1, 1],
                "facility": [0, 1, 0, 1],
                "cost": [1.0, 100.0, 100.0, 1.0],
            }
        )

        result = solve_facility_location(
            customer_data, facility_data, transportation_cost
        )

        # Check that result has expected structure
        self.assertIn("solution_value", result)
        self.assertIn("facilities_opened", result)
        self.assertIn("assignments", result)

        # Check solution value is close to expected
        self.assertAlmostEqual(result["solution_value"], 30.0, places=5)

        # Check both facilities are opened
        opened = result["facilities_opened"]
        self.assertIsInstance(opened, pd.Series)
        self.assertEqual(opened.sum(), 2)

        # Check assignments
        assignments = result["assignments"]
        self.assertIsInstance(assignments, pd.DataFrame)
        self.assertTrue(
            all(
                col in assignments.columns
                for col in ["customer", "facility", "assignment"]
            )
        )

    def test_all_demand_satisfied(self):
        """Test that all customer demand is satisfied."""
        customer_data = pd.DataFrame(
            {"customer": [0, 1, 2], "demand": [10.0, 15.0, 20.0]}
        )

        facility_data = pd.DataFrame(
            {"facility": [0, 1], "capacity": [30.0, 30.0], "fixed_cost": [50.0, 50.0]}
        )

        transportation_cost = pd.DataFrame(
            {
                "customer": [0, 0, 1, 1, 2, 2],
                "facility": [0, 1, 0, 1, 0, 1],
                "cost": [1.0, 2.0, 2.0, 1.0, 3.0, 1.0],
            }
        )

        result = solve_facility_location(
            customer_data, facility_data, transportation_cost
        )

        # Sum of assignments for each customer should equal their demand
        assignments = result["assignments"]
        total_assigned = assignments.groupby("customer")["assignment"].sum()

        for customer in customer_data["customer"]:
            self.assertAlmostEqual(
                total_assigned.loc[customer],
                customer_data.loc[customer_data["customer"] == customer, "demand"].iloc[
                    0
                ],
                places=5,
            )

    def test_capacity_constraints_respected(self):
        """Test that facility capacity constraints are respected."""
        customer_data = pd.DataFrame(
            {"customer": [0, 1, 2], "demand": [10.0, 10.0, 10.0]}
        )

        facility_data = pd.DataFrame(
            {"facility": [0, 1], "capacity": [15.0, 20.0], "fixed_cost": [10.0, 20.0]}
        )

        transportation_cost = pd.DataFrame(
            {
                "customer": [0, 0, 1, 1, 2, 2],
                "facility": [0, 1, 0, 1, 0, 1],
                "cost": [1.0, 2.0, 1.0, 2.0, 1.0, 2.0],
            }
        )

        result = solve_facility_location(
            customer_data, facility_data, transportation_cost
        )

        # Sum of assignments to each facility should not exceed capacity
        assignments = result["assignments"]
        total_at_facility = assignments.groupby("facility")["assignment"].sum()

        for facility in facility_data["facility"]:
            capacity = facility_data.loc[
                facility_data["facility"] == facility, "capacity"
            ].iloc[0]
            if facility in total_at_facility.index:
                self.assertLessEqual(total_at_facility.loc[facility], capacity + 1e-5)

    def test_assignment_only_to_open_facilities(self):
        """Test that customers are only assigned to open facilities."""
        customer_data = pd.DataFrame({"customer": [0, 1], "demand": [5.0, 5.0]})

        facility_data = pd.DataFrame(
            {
                "facility": [0, 1, 2],
                "capacity": [10.0, 10.0, 10.0],
                "fixed_cost": [5.0, 100.0, 100.0],  # Facility 0 much cheaper
            }
        )

        transportation_cost = pd.DataFrame(
            {
                "customer": [0, 0, 0, 1, 1, 1],
                "facility": [0, 1, 2, 0, 1, 2],
                "cost": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            }
        )

        result = solve_facility_location(
            customer_data, facility_data, transportation_cost
        )

        # Get facilities with non-zero assignments
        assignments = result["assignments"]
        assigned_facilities = assignments[assignments["assignment"] > 1e-6][
            "facility"
        ].unique()

        # Check all assigned facilities are open
        opened = result["facilities_opened"]
        for facility in assigned_facilities:
            self.assertTrue(opened.loc[facility] > 0.5)  # Binary variable


class TestFacilityLocationEdgeCases(unittest.TestCase):
    """Test edge cases and error handling."""

    def test_single_customer_single_facility(self):
        """Test trivial case with one customer and one facility."""
        customer_data = pd.DataFrame({"customer": [0], "demand": [10.0]})

        facility_data = pd.DataFrame(
            {"facility": [0], "capacity": [10.0], "fixed_cost": [5.0]}
        )

        transportation_cost = pd.DataFrame(
            {"customer": [0], "facility": [0], "cost": [2.0]}
        )

        result = solve_facility_location(
            customer_data, facility_data, transportation_cost
        )

        # Facility must be opened
        self.assertAlmostEqual(result["facilities_opened"].iloc[0], 1.0)
        # Total cost = fixed + transportation
        self.assertAlmostEqual(result["solution_value"], 5.0 + 10.0 * 2.0)

    def test_zero_fixed_costs(self):
        """Test when all facilities have zero fixed cost."""
        customer_data = pd.DataFrame({"customer": [0, 1], "demand": [5.0, 5.0]})

        facility_data = pd.DataFrame(
            {"facility": [0, 1], "capacity": [10.0, 10.0], "fixed_cost": [0.0, 0.0]}
        )

        transportation_cost = pd.DataFrame(
            {
                "customer": [0, 0, 1, 1],
                "facility": [0, 1, 0, 1],
                "cost": [1.0, 2.0, 2.0, 1.0],
            }
        )

        result = solve_facility_location(
            customer_data, facility_data, transportation_cost
        )

        # Should succeed and return valid solution
        self.assertIsNotNone(result)
        self.assertIn("solution_value", result)

    def test_insufficient_capacity_infeasible(self):
        """Test that infeasible problem is detected."""
        customer_data = pd.DataFrame({"customer": [0], "demand": [100.0]})

        facility_data = pd.DataFrame(
            {
                "facility": [0, 1],
                "capacity": [10.0, 20.0],  # Total capacity = 30 < demand = 100
                "fixed_cost": [5.0, 5.0],
            }
        )

        transportation_cost = pd.DataFrame(
            {"customer": [0, 0], "facility": [0, 1], "cost": [1.0, 1.0]}
        )

        with self.assertRaises(ValueError) as context:
            solve_facility_location(customer_data, facility_data, transportation_cost)

        self.assertIn("infeasible", str(context.exception).lower())


class TestFacilityLocationDataValidation(unittest.TestCase):
    """Test input validation and error handling."""

    def test_missing_customer_columns(self):
        """Test error when customer data is missing required columns."""
        customer_data = pd.DataFrame(
            {
                "customer": [0, 1],
                # Missing 'demand' column
            }
        )

        facility_data = pd.DataFrame(
            {"facility": [0, 1], "capacity": [10.0, 10.0], "fixed_cost": [5.0, 5.0]}
        )

        transportation_cost = pd.DataFrame(
            {
                "customer": [0, 0, 1, 1],
                "facility": [0, 1, 0, 1],
                "cost": [1.0, 1.0, 1.0, 1.0],
            }
        )

        with self.assertRaises(ValueError):
            solve_facility_location(customer_data, facility_data, transportation_cost)

    def test_missing_facility_columns(self):
        """Test error when facility data is missing required columns."""
        customer_data = pd.DataFrame({"customer": [0, 1], "demand": [10.0, 10.0]})

        facility_data = pd.DataFrame(
            {
                "facility": [0, 1],
                # Missing 'capacity' and 'fixed_cost' columns
            }
        )

        transportation_cost = pd.DataFrame(
            {
                "customer": [0, 0, 1, 1],
                "facility": [0, 1, 0, 1],
                "cost": [1.0, 1.0, 1.0, 1.0],
            }
        )

        with self.assertRaises(ValueError):
            solve_facility_location(customer_data, facility_data, transportation_cost)

    def test_missing_transportation_columns(self):
        """Test error when transportation cost data is missing required columns."""
        customer_data = pd.DataFrame({"customer": [0, 1], "demand": [10.0, 10.0]})

        facility_data = pd.DataFrame(
            {"facility": [0, 1], "capacity": [10.0, 10.0], "fixed_cost": [5.0, 5.0]}
        )

        transportation_cost = pd.DataFrame(
            {
                "customer": [0, 0, 1, 1],
                "facility": [0, 1, 0, 1],
                # Missing 'cost' column
            }
        )

        with self.assertRaises(ValueError):
            solve_facility_location(customer_data, facility_data, transportation_cost)

    def test_negative_demand(self):
        """Test error when demand is negative."""
        customer_data = pd.DataFrame(
            {"customer": [0, 1], "demand": [-10.0, 10.0]}  # Negative demand
        )

        facility_data = pd.DataFrame(
            {"facility": [0, 1], "capacity": [10.0, 10.0], "fixed_cost": [5.0, 5.0]}
        )

        transportation_cost = pd.DataFrame(
            {
                "customer": [0, 0, 1, 1],
                "facility": [0, 1, 0, 1],
                "cost": [1.0, 1.0, 1.0, 1.0],
            }
        )

        with self.assertRaises(ValueError):
            solve_facility_location(customer_data, facility_data, transportation_cost)

    def test_negative_capacity(self):
        """Test error when capacity is negative."""
        customer_data = pd.DataFrame({"customer": [0, 1], "demand": [10.0, 10.0]})

        facility_data = pd.DataFrame(
            {
                "facility": [0, 1],
                "capacity": [-10.0, 10.0],  # Negative capacity
                "fixed_cost": [5.0, 5.0],
            }
        )

        transportation_cost = pd.DataFrame(
            {
                "customer": [0, 0, 1, 1],
                "facility": [0, 1, 0, 1],
                "cost": [1.0, 1.0, 1.0, 1.0],
            }
        )

        with self.assertRaises(ValueError):
            solve_facility_location(customer_data, facility_data, transportation_cost)

    def test_negative_transportation_cost(self):
        """Test error when transportation cost is negative."""
        customer_data = pd.DataFrame({"customer": [0, 1], "demand": [10.0, 10.0]})

        facility_data = pd.DataFrame(
            {"facility": [0, 1], "capacity": [10.0, 10.0], "fixed_cost": [5.0, 5.0]}
        )

        transportation_cost = pd.DataFrame(
            {
                "customer": [0, 0, 1, 1],
                "facility": [0, 1, 0, 1],
                "cost": [-1.0, 1.0, 1.0, 1.0],  # Negative cost
            }
        )

        with self.assertRaises(ValueError):
            solve_facility_location(customer_data, facility_data, transportation_cost)

    def test_empty_dataframes(self):
        """Test error when input dataframes are empty."""
        customer_data = pd.DataFrame({"customer": [], "demand": []})

        facility_data = pd.DataFrame({"facility": [], "capacity": [], "fixed_cost": []})

        transportation_cost = pd.DataFrame({"customer": [], "facility": [], "cost": []})

        with self.assertRaises(ValueError):
            solve_facility_location(customer_data, facility_data, transportation_cost)


if __name__ == "__main__":
    unittest.main()
