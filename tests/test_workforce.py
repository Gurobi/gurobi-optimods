import unittest

from pandas.api.types import is_object_dtype, is_datetime64_any_dtype, is_numeric_dtype

from gurobi_optimods.workforce import solve_workforce_scheduling
from gurobi_optimods.datasets import load_workforce


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
