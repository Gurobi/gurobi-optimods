"""
Diet problem with a pandas interface
"""

from dataclasses import dataclass

import gurobipy as gp
from gurobipy import GRB
import gurobipy_pandas as gppd
import pandas as pd


@dataclass
class DietResult:
    menu: pd.Series
    total_cost: float


def solve_diet_problem(categories, foods, values):
    """
    A sphinx-compatible docstring

    :param data1: Data structure for first argument
    :type data1: pd.DataFrame
    """
    with gp.Env() as env, gp.Model(env=env) as model:
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
