"""
Diet Problem
------------
"""

from dataclasses import dataclass

import gurobipy as gp
from gurobipy import GRB
import gurobipy_pandas as gppd
import pandas as pd

from gurobi_optimods.utils import optimod


@dataclass
class DietResult:
    menu: pd.Series
    total_cost: float


@optimod()
def solve_diet_problem(categories, foods, values, *, create_env):
    """
    Choose quantities of foods to eat in order to meet the required
    nutrient amounts in each category, for minimum total cost.

    Parameters
    ----------
    categories : DataFrame
        Dataframe with columns (category, min, max)
    foods : DataFrame
        Dataframe with columns (food, cost)
    values : DataFrame
        Dataframe with columns (food, category, value)

    Returns
    -------
    DataFrame
        A menu!

    Raises
    ------
    ValueError
        If a feasible diet isn't possible

    """
    with create_env() as env, gp.Model(env=env) as model:
        # Build the model
        quantity = gppd.add_vars(
            model, foods.set_index("food"), obj="cost", name="quantity"
        )
        amounts = (
            values.join(quantity, on="food")
            .assign(amount=lambda df: df["value"] * df["quantity"])
            .groupby("category")["amount"]
            .sum()
        )
        (
            categories.join(amounts, on="category")
            .gppd.add_constrs(model, "amount >= min", name="lower")
            .gppd.add_constrs(model, "amount <= max", name="upper")
        )
        # Solve, post-process and return solution
        model.optimize()
        if model.Status == GRB.INFEASIBLE:
            raise ValueError("Unsatisfiable diet")
        return DietResult(menu=quantity.gppd.X, total_cost=model.ObjVal)
