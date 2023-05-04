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
            set(data.keys()), {"preferences", "shift_requirements", "worker_limits"}
        )

        self.assertEqual(
            set(data.preferences.columns), {"Worker", "Shift", "Preference"}
        )
        self.assertTrue(is_object_dtype(data.preferences["Worker"]))
        self.assertTrue(is_numeric_dtype(data.preferences["Preference"]))
        self.assertTrue(is_datetime64_any_dtype(data.preferences["Shift"]))

        self.assertEqual(set(data.shift_requirements.columns), {"Shift", "Required"})
        self.assertTrue(is_datetime64_any_dtype(data.shift_requirements["Shift"]))
        self.assertTrue(is_numeric_dtype(data.shift_requirements["Required"]))

        self.assertEqual(
            set(data.worker_limits.columns),
            {"Worker", "MinShifts", "MaxShifts"},
        )
        self.assertTrue(is_numeric_dtype(data.worker_limits["MinShifts"]))
        self.assertTrue(is_numeric_dtype(data.worker_limits["MaxShifts"]))

    def test_no_option(self):
        # Simple example where there is only one way to cover requirements
        preferences = read_csv(
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
            preferences=preferences, shift_requirements=shift_requirements
        )

        self.assertIsInstance(assignments, pd.DataFrame)
        self.assertIsNot(assignments, preferences)
        assert_frame_equal(assignments, preferences)

    def test_preferences(self):
        # Choose an assignment which maximises sum of preferences
        preferences = read_csv(
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
            preferences=preferences, shift_requirements=shift_requirements
        )

        expected = read_csv(
            """
            Worker,Shift,Preference
            Alice,2022-07-03,2.0
            Bob,2022-07-02,2.0
            """
        )
        self.assertIsInstance(assignments, pd.DataFrame)
        self.assertIsNot(assignments, preferences)
        assert_frame_equal(
            assignments.sort_values(["Worker", "Shift"]).reset_index(drop=True),
            expected,
        )

    def test_constrained(self):
        # Test min/max shift constraints
        preferences = read_csv(
            """
            Worker,Shift,Preference
            Alice,2022-07-01,1.0
            Alice,2022-07-02,1.0
            Alice,2022-07-03,1.0
            Bob,2022-07-01,2.0
            Bob,2022-07-02,2.0
            Bob,2022-07-03,2.0
            Joy,2022-07-01,3.0
            Joy,2022-07-02,3.1
            Joy,2022-07-03,3.1
            """
        )
        shift_requirements = read_csv(
            """
            Shift,Required
            2022-07-01,1
            2022-07-02,1
            2022-07-03,1
            """
        )

        with self.subTest("upperlimits"):

            worker_limits = read_csv(
                """
                Worker,MinShifts,MaxShifts
                Alice,0,3
                Bob,0,1
                Joy,0,2
                """
            )

            assignments = solve_workforce_scheduling(
                preferences=preferences,
                shift_requirements=shift_requirements,
                worker_limits=worker_limits,
            )
            expected = read_csv(
                """
                Worker,Shift,Preference
                Bob,2022-07-01,2.0
                Joy,2022-07-02,3.1
                Joy,2022-07-03,3.1
                """
            )
            self.assertIsInstance(assignments, pd.DataFrame)
            self.assertIsNot(assignments, preferences)
            assert_frame_equal(
                assignments.sort_values(["Shift"]).reset_index(drop=True),
                expected,
            )

        with self.subTest("lowerlimits"):

            worker_limits = read_csv(
                """
                Worker,MinShifts,MaxShifts
                Alice,1,3
                Bob,0,3
                Joy,0,3
                """
            )

            assignments = solve_workforce_scheduling(
                preferences=preferences,
                shift_requirements=shift_requirements,
                worker_limits=worker_limits,
            )
            expected = read_csv(
                """
                Worker,Shift,Preference
                Alice,2022-07-01,1.0
                Joy,2022-07-02,3.1
                Joy,2022-07-03,3.1
                """
            )
            self.assertIsInstance(assignments, pd.DataFrame)
            self.assertIsNot(assignments, preferences)
            assert_frame_equal(
                assignments.sort_values(["Shift"]).reset_index(drop=True),
                expected,
            )

    def test_rolling(self):

        preferences = read_csv(
            """
            Worker,Shift,Preference
            Alice,2022-07-01,1.0
            Alice,2022-07-02,2.0
            Alice,2022-07-03,3.0
            Alice,2022-07-04,4.0
            Bob,2022-07-01,5.1
            Bob,2022-07-02,6.2
            Bob,2022-07-03,7.3
            Bob,2022-07-04,8.4
            """
        ).assign(Shift=lambda df: pd.to_datetime(df["Shift"]))
        shift_requirements = read_csv(
            """
            Shift,Required
            2022-07-01,1
            2022-07-02,1
            2022-07-03,1
            2022-07-04,1
            """
        ).assign(Shift=lambda df: pd.to_datetime(df["Shift"]))

        assignments = solve_workforce_scheduling(
            preferences=preferences,
            shift_requirements=shift_requirements,
            rolling_window=pd.Timedelta(days=2),
            rolling_limit=1,
        )

        expected = read_csv(
            """
            Worker,Shift,Preference
            Alice,2022-07-01,1.0
            Bob,2022-07-02,6.2
            Alice,2022-07-03,3.0
            Bob,2022-07-04,8.4
            """
        ).assign(Shift=lambda df: pd.to_datetime(df["Shift"]))
        self.assertIsInstance(assignments, pd.DataFrame)
        self.assertIsNot(assignments, preferences)
        assert_frame_equal(
            assignments.sort_values(["Shift"]).reset_index(drop=True),
            expected,
        )

    def test_infeasibility(self):
        # Should raise an exception
        data = load_workforce()

        with self.assertRaises(ValueError):
            solve_workforce_scheduling(
                data.preferences,
                data.shift_requirements,
                data.worker_limits,
                rolling_window=pd.Timedelta("4D"),
                rolling_limit=3,
            )
