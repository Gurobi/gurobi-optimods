import gurobipy as gp
from gurobipy import GRB
import numpy as np
import math

from gurobi_optimods.utils import optimod


@optimod()
def max_sharpe_ratio(Q: np.ndarray, mu: np.ndarray, rf_rate: float, *, create_env):
    """
    Solve the problem of finding a portfolio that maximizes the
    Sharpe ratio.

    :param Q: 2D positive-semidefinite variance-covariance matrix.
    :type Q: :class:`np.ndarray `
    :param mu: Return rates.
    :type mu: :class:`np.ndarray`
    :param rf_rate: Risk-free rate of return.
    :type rf_rate: float
    :return: Portfolio that maximizes the Sharpe ratio.
    :rtype: :class:`np.ndarray`
    :param silent: Optional. Boolean with whether output should be printed.
    :type silent: :class:`bool`
    :param logfile: Optional. String with file path with logger and Gurobi.
    :type logfile: :class:`str`
    :return: Sharpe ratio of the portfolio.
    :rtype: float
    """
    with create_env() as env, gp.Model("sharpe_ratio", env=env) as model:
        y = model.addMVar(mu.size, name="y")
        model.addConstr((mu - rf_rate) @ y == 1)
        model.setObjective(y @ Q @ y, sense=GRB.MINIMIZE)

        model.optimize()

        # Translate solution to original variable space before returning
        return y.X / y.X.sum(), 1 / math.sqrt(model.ObjVal)
