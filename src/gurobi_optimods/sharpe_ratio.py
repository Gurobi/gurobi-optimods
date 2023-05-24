import gurobipy as gp
from gurobipy import GRB
import numpy as np
import pandas as pd
import math

from gurobi_optimods.utils import optimod


@optimod()
def max_sharpe_ratio(Q, mu, rf_rate=0, *, create_env):
    """
    Solve the problem of finding a portfolio that maximizes the
    Sharpe ratio.

    :param Q: 2D positive-semidefinite variance-covariance matrix.
    :type Q: :class:`np.ndarray|pd.DataFrame`
    :param mu: Return rates.
    :type mu: :class:`np.ndarray|pd.Series`
    :param rf_rate: Risk-free rate of return (optional, defaults to ``0``).
    :type rf_rate: float
    :return: Portfolio that maximizes the Sharpe ratio.
    :rtype: :class:`np.ndarray|pd.Series`
    :param silent: Optional. Boolean with whether output should be printed.
    :type silent: :class:`bool`
    :param logfile: Optional. String with file path with logger and Gurobi.
    :type logfile: :class:`str`
    :return: Optimal portfolio.
    :rtype: :class:`np.ndarray|pd.Series`
    :return: Sharpe ratio of the portfolio.
    :rtype: float
    """
    indices = None

    if isinstance(Q, pd.DataFrame):
        indices = Q.index
        Q = Q.to_numpy()
    elif not isinstance(Q, np.ndarray):
        raise ValueError(f"Unknown covariance matrix type: {type(Q)}")

    if isinstance(mu, pd.Series):
        if indices is None:
            indices = mu.index
        elif not mu.index.equals(indices):
            raise ValueError("Indices of Q and mu are misaligned")
        mu = mu.to_numpy()
    elif not isinstance(mu, np.ndarray):
        raise ValueError(f"Unknown return rates type: {type(mu)}")

    if not isinstance(rf_rate, float) and not isinstance(rf_rate, int):
        raise ValueError(f"Unknown risk-free return rate type: {type(rf_rate)}")
    elif (mu < rf_rate).all():
        raise ValueError(
            f"No expected returns are greater than risk-free return rate of {rf_rate}"
        )

    portfolio, ratio = _max_sharpe_ratio_numpy(Q, mu, rf_rate, create_env)

    if indices is not None:
        portfolio = pd.Series(data=portfolio, index=indices)

    return portfolio, ratio


def _max_sharpe_ratio_numpy(Q, mu, rf_rate, create_env):
    with create_env() as env, gp.Model("sharpe_ratio", env=env) as model:
        y = model.addMVar(mu.size, name="y")
        model.addConstr((mu - rf_rate) @ y == 1)
        model.setObjective(y @ Q @ y, sense=GRB.MINIMIZE)

        model.optimize()

        # Translate solution to original variable space before returning
        return y.X / y.X.sum(), 1 / math.sqrt(model.ObjVal)
