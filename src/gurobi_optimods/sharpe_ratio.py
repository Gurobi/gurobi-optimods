"""
Sharpe Ratio
------------
"""

import gurobipy as gp
from gurobipy import GRB
import numpy as np
import pandas as pd
import math

from gurobi_optimods.utils import optimod


@optimod()
def max_sharpe_ratio(cov_matrix, mu, rf_rate=0, *, create_env):
    """
    Solve the problem of finding a portfolio that maximizes the
    Sharpe ratio.

    :param cov_matrix: 2D positive-semidefinite variance-covariance matrix :math:`\Sigma`
    :type cov_matrix: :class:`np.ndarray|pd.DataFrame`
    :param mu: Expected return rates :math:`\mu`
    :type mu: :class:`np.ndarray|pd.Series`
    :param rf_rate: Non-negative risk-free rate of return (optional, defaults to ``0``)
    :type rf_rate: :class:`float` >= 0
    :return: Portfolio that maximizes the Sharpe ratio
    :rtype: :class:`np.ndarray|pd.Series`
    :return: Sharpe ratio of the portfolio
    :rtype: :class:`float`
    """
    indices = None

    if isinstance(cov_matrix, pd.DataFrame):
        indices = cov_matrix.index
        cov_matrix = cov_matrix.to_numpy()
    elif not isinstance(cov_matrix, np.ndarray):
        raise TypeError(f"Unknown covariance matrix type: {type(cov_matrix)}")

    if cov_matrix.ndim != 2:
        raise ValueError(
            f"Covariance matrix should be in 2 dimensions, not {cov_matrix.ndim}"
        )

    if isinstance(mu, pd.Series):
        if indices is None:
            indices = mu.index
        elif not mu.index.equals(indices):
            raise ValueError("Indices of cov_matrix and mu are misaligned")
        mu = mu.to_numpy()
    elif not isinstance(mu, np.ndarray):
        raise TypeError(f"Unknown return rates type: {type(mu)}")

    if mu.ndim != 1:
        raise ValueError(f"Return rates should be in 1 dimension, not {mu.ndim}")

    if not isinstance(rf_rate, float) and not isinstance(rf_rate, int):
        raise TypeError(f"Unknown risk-free return rate type: {type(rf_rate)}")
    elif rf_rate < 0:
        raise ValueError("Risk-free return rate should be non-negative")
    elif (mu < rf_rate).all():
        raise ValueError(
            f"No expected returns are greater than risk-free return rate of {rf_rate}"
        )

    portfolio, ratio = _max_sharpe_ratio_numpy(cov_matrix, mu, rf_rate, create_env)

    if indices is not None:
        portfolio = pd.Series(data=portfolio, index=indices)

    return portfolio, ratio


def _max_sharpe_ratio_numpy(cov_matrix, mu, rf_rate, create_env):
    with create_env() as env, gp.Model("sharpe_ratio", env=env) as model:
        y = model.addMVar(mu.size, name="y")
        model.addConstr((mu - rf_rate) @ y == 1)
        model.setObjective(y @ cov_matrix @ y, sense=GRB.MINIMIZE)

        model.optimize()

        # Translate solution to original variable space before returning
        return y.X / y.X.sum(), 1 / math.sqrt(model.ObjVal)
