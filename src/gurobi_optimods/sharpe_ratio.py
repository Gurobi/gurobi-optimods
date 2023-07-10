"""
Sharpe Ratio
------------
"""

import math
from dataclasses import dataclass
from typing import Union

import gurobipy as gp
import numpy as np
import pandas as pd
from gurobipy import GRB

from gurobi_optimods.utils import optimod


@optimod()
def max_sharpe_ratio(cov_matrix, mu, rf_rate=0, *, create_env):
    """
    Solve the problem of finding a portfolio that maximizes the
    Sharpe ratio.

    Parameters
    ----------
    cov_matrix : ndarray or DataFrame
        2D positive-semidefinite variance-covariance matrix :math:`\Sigma`
    mu : ndarray or Series
        Expected return rates :math:`\mu`
    rf_rate : float >= 0, optional
        Non-negative risk-free rate of return (defaults to ``0``)

    Returns
    -------
    result : SharpeRatioResult
        A data class representing the portfolio that maximizes the Sharpe
        ratio:

        * ``result.x``: The relative portfolio allocations :math:`x`
        * ``result.sharpe_ratio``: The Sharpe ratio of the portfolio
        * ``result.ret``: The estimated return :math:`\mu^T x` of the portfolio
        * ``result.risk``: The estimated risk :math:`x^T \Sigma x` of the
          portfolio
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

    result = _max_sharpe_ratio_numpy(cov_matrix, mu, rf_rate, create_env)

    if indices is not None:
        result.x = pd.Series(data=result.x, index=indices)

    return result


def _max_sharpe_ratio_numpy(cov_matrix, mu, rf_rate, create_env):
    with create_env() as env, gp.Model("sharpe_ratio", env=env) as model:
        y = model.addMVar(mu.size, name="y")
        model.addConstr((mu - rf_rate) @ y == 1)
        model.setObjective(y @ cov_matrix @ y, sense=GRB.MINIMIZE)

        model.optimize()

        # Translate solution to original variable space
        x = y.X / y.X.sum()
        ret = mu @ x
        risk = x @ cov_matrix @ x
        sharpe_ratio = (ret - rf_rate) / math.sqrt(risk)
        return SharpeRatioResult(x, sharpe_ratio, ret, risk)


@dataclass
class SharpeRatioResult:
    """
    Data class representing the portfolio that maximizes the Sharpe ratio.


    Attributes
    ----------
    x : ndarray or Series
        The relative portfolio allocations :math:`x`
    sharpe_ratio : float
        The Sharpe ratio of the portfolio
    ret : float
        The estimated return :math:`\mu^T x` of the portfolio
    risk : float
        The estimated risk :math:`x^T \Sigma x` of the portfolio
    """

    x: Union[np.ndarray, pd.Series]
    sharpe_ratio: float
    ret: float
    risk: float
