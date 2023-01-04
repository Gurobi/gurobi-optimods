# Unit tests for your mod. This is separate from the example tests, it should
# test smaller pieces of functionality, and can also test multiple cases.
# This template should be copied to tests/test_<mod-name>.py

import unittest

import pandas as pd

from gurobi_optimods.datasets import load_diet
from gurobi_optimods.diet import solve_diet_problem


class TestDiet(unittest.TestCase):
    def test_datasets(self):
        # If you added something to optimods.datasets, test it here
        data = load_diet()
        self.assertEqual(set(data.categories.columns), {"category", "min", "max"})
        self.assertEqual(set(data.foods.columns), {"food", "cost"})
        self.assertEqual(
            set(data.nutrition_values.columns), {"food", "category", "value"}
        )

    def test_survive_on_hamburgers(self):
        # I have one option. How many hamburgers must I eat to survive?
        categories = pd.DataFrame(
            {
                "category": ["calories"],
                "min": [1800],
                "max": [2200],
            }
        )
        foods = pd.DataFrame({"food": ["hamburger"], "cost": [2.5]})
        values = pd.DataFrame(
            {
                "food": ["hamburger"],
                "category": ["calories"],
                "value": [400],
            }
        )

        # User passes the dataframes, gets a result.
        diet = solve_diet_problem(categories=categories, foods=foods, values=values)

        # Expect the result is an object containing a menu (series with required
        # quantities) and total cost.
        self.assertIsInstance(diet.menu, pd.Series)
        self.assertEqual(set(diet.menu.index), {"hamburger"})
        self.assertEqual(diet.menu["hamburger"], 1800 / 400)
        self.assertEqual(diet.total_cost, 1800 / 400 * 2.5)

    def test_kryptonite_hamburgers(self):
        # Can I eat enough hamburgers to survive without losing my superpowers?
        categories = pd.DataFrame(
            {
                "category": ["calories", "kryptonite"],
                "min": [1800, 0],
                "max": [2200, 10],
            }
        )
        foods = pd.DataFrame({"food": ["hamburger"], "cost": [2.5]})
        values = pd.DataFrame(
            {
                "food": ["hamburger", "hamburger"],
                "category": ["calories", "kryptonite"],
                "value": [400, 4],
            }
        )

        # Eating enough calories exceeds my kryptonite tolerance: infeasible.
        with self.assertRaisesRegex(ValueError, "Unsatisfiable diet"):
            diet = solve_diet_problem(categories=categories, foods=foods, values=values)

    def test_healthy_alternative(self):
        # There is an expensive kryptonite-free hamburger alternative on the market.
        categories = pd.DataFrame(
            {
                "category": ["calories", "kryptonite"],
                "min": [1800, 0],
                "max": [2200, 10],
            }
        )
        foods = pd.DataFrame({"food": ["hamburger", "broccoli"], "cost": [2.5, 4.0]})
        values = pd.DataFrame(
            {
                "food": ["hamburger", "hamburger", "broccoli", "broccoli"],
                "category": ["calories", "kryptonite", "calories", "kryptonite"],
                "value": [400, 4, 200, 0],
            }
        )

        # User passes the dataframes, gets a result.
        diet = solve_diet_problem(categories=categories, foods=foods, values=values)

        # Use the cheap calories first, up to my kryptonite limit.
        self.assertIsInstance(diet.menu, pd.Series)
        self.assertEqual(set(diet.menu.index), {"hamburger", "broccoli"})
        self.assertEqual(diet.menu["hamburger"], 10 / 4)
        self.assertEqual(diet.menu["broccoli"], (1800 - (10 / 4) * 400) / 200)
        self.assertEqual(diet.total_cost, 4.0 * 4.0 + 2.5 * 2.5)
