import unittest
from io import StringIO
from textwrap import dedent

import pandas as pd
from pandas.api.types import is_object_dtype, is_datetime64_any_dtype, is_numeric_dtype
from pandas.testing import assert_frame_equal

from gurobi_optimods.workforce import solve_workforce_scheduling
from gurobi_optimods.datasets import load_workforce


def read_csv(text):
    return pd.read_csv(StringIO(dedent(text)))


class TestWorkforceScheduling(unittest.TestCase):
    def test_dataset(self):
        data = load_workforce()
        self.assertEqual(
            set(data.keys()), {"availability", "pay_rates", "shift_requirements"}
        )

        self.assertEqual(set(data.availability.columns), {"Worker", "Shift"})
        self.assertTrue(is_object_dtype(data.availability["Worker"]))
        self.assertTrue(is_datetime64_any_dtype(data.availability["Shift"]))

        self.assertEqual(set(data.pay_rates.columns), {"Worker", "PayRate"})
        self.assertTrue(is_object_dtype(data.pay_rates["Worker"]))
        self.assertTrue(is_numeric_dtype(data.pay_rates["PayRate"]))

        self.assertEqual(set(data.shift_requirements.columns), {"Shift", "Required"})
        self.assertTrue(is_datetime64_any_dtype(data.shift_requirements["Shift"]))
        self.assertTrue(is_numeric_dtype(data.shift_requirements["Required"]))

    def test_no_option(self):
        # Simple example where there is only one way to cover requirements
        availability = read_csv(
            """
            Worker,Shift,Preference
            Bob,2022-07-02,1.0
            Alice,2022-07-03,1.0
            """
        )
        shift_requirements = read_csv(
            """
            Shift,Required
            2022-07-02,1
            2022-07-03,1
            """
        )

        assignments = solve_workforce_scheduling(
            availability=availability, shift_requirements=shift_requirements
        )

        self.assertIsInstance(assignments, pd.DataFrame)
        self.assertIsNot(assignments, availability)
        assert_frame_equal(assignments, availability)

    def test_preferences(self):
        # Choose an assignment which maximises preferences
        availability = read_csv(
            """
            Worker,Shift,Preference
            Alice,2022-07-02,1.0
            Alice,2022-07-03,2.0
            Bob,2022-07-02,2.0
            Bob,2022-07-03,1.0
            """
        )
        shift_requirements = read_csv(
            """
            Shift,Required
            2022-07-02,1
            2022-07-03,1
            """
        )

        assignments = solve_workforce_scheduling(
            availability=availability, shift_requirements=shift_requirements
        )

        expected = read_csv(
            """
            Worker,Shift,Preference
            Alice,2022-07-03,2.0
            Bob,2022-07-02,2.0
            """
        )
        self.assertIsInstance(assignments, pd.DataFrame)
        self.assertIsNot(assignments, availability)
        assert_frame_equal(
            assignments.sort_values(["Worker", "Shift"]).reset_index(drop=True),
            expected,
        )
