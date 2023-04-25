import gurobipy as gp
from gurobipy import GRB
import numpy as np
import scipy.sparse as sp
import pandas as pd


class Portfolio:
    def __init__(
        self,
        H,
        Sigma,
        mu,
        initial_holdings=None,
        maxnum_transactions=0,
        maxnum_assets=0,
        min_invest_long=0,
        min_invest_short=0,
        leverage=0,
    ) -> None:

        self.factor_matrix = H  # H.T@H is the covariance matrix
        self.covariance = Sigma  # Sigma is the covariance matrix
        self.mu = mu  # estimated first moments of return function
        self.initial_holdings = initial_holdings  # existing allocation
        self.maxnum_transactions = maxnum_transactions  # No more than 20 trades
        self.maxnum_assets = (
            maxnum_assets  # No more than 50 assets at a time (e.g., online problem)
        )
        self.min_invest_long = min_invest_long  # For long allocations, need at least 5% of total investment
        self.min_invest_short = min_invest_short  # For short allocations, need at least 1% of total investment
        self.leverage = leverage  # Maximum 10% short selling


class MeanVariancePortfolio:
    """Optimal mean-variance portfolio solver.

    Instantiate an object of :class:`MeanVariancePortfolio` for given
    covariance matrix and return vector.  Use
    :meth:`MeanVariancePortfolio.efficient_portfolio` to solve for efficient
    portfolios with given parameters.

    :param Sigma: Covariance matrix
    :type Sigma: 2-d :class:`np.ndarray`
    :param mu: Return vector
    :type mu: 1-d :class:`np.ndarray`

    """

    def __init__(
        self,
        Sigma,
        mu,
    ) -> None:
        if isinstance(Sigma, pd.DataFrame):
            self.resultType = "pandas"
            self.index = Sigma.index
            self.covariance = Sigma.to_numpy()
        elif isinstance(Sigma, np.ndarray):
            self.covariance = Sigma
            self.resultType = "numpy"
        else:
            raise TypeError("Incompatible type of Sigma")

        if isinstance(mu, pd.Series):
            self.resultType = "pandas"
            self.mu = mu.to_numpy()
        elif isinstance(mu, np.ndarray):
            self.mu = mu
            self.resultType = "numpy"
        else:
            raise TypeError("Incompatible type of mu")

    def _convert_result(self, x):
        if self.resultType == "numpy":
            return x
        elif self.resultType == "pandas":
            return pd.Series(x, index=self.index)
        else:
            assert False

    def minimize_risk(self, expected_return):
        with gp.Env() as env, gp.Model("min_risk", env=env) as m:
            x = m.addMVar(shape=self.mu.shape, name="x")

            m.addConstr(x.sum() == 1, name="fully_invested")
            m.addConstr(self.mu @ x >= expected_return, name="expected_return")
            m.setObjective(x @ self.covariance @ x)

            m.optimize()

            if m.Status == GRB.OPTIMAL:
                return self._convert_result(x.X)

    def maximize_return(self, max_risk):
        with gp.Env() as env, gp.Model("max_return", env=env) as m:
            x = m.addMVar(shape=self.mu.shape, name="x")
            m.addConstr(x.sum() == 1, name="fully_invested")
            m.addConstr(x @ self.covariance @ x <= max_risk, name="max_risk")
            m.setObjective(self.mu @ x, GRB.MAXIMIZE)

            m.optimize()

            if m.Status == GRB.OPTIMAL:
                return self._convert_result(x.X)

    def efficient_portfolio(
        self, gamma, max_trades=None, fees_buy=None, min_buy_in=None
    ):
        """
        Compute efficient portfolio for given paramters

        :param gamma: Risk aversion cofficient for balancing risk and return;
            the resulting objective functions is
            :math:`\mu^T x - 0.5 \gamma x^T \Sigma x`
        :type gamma: :class:`float` >= 0
        :param max_trades: Upper limit on the number of trades
        :type max_trades: :class:`int` >= 0
        :param fees_buys: Fixed-charge cost for each traded position, relative
            to total portfolio value
        :type fees_buys: :class:`float` >= 0
        :param min_buy_in: Lower bound on the volume on a traded long position,
            relative to total portfolio value
        :type min_buy_in: :class:`float` >= 0
        """
        with gp.Env() as env, gp.Model("efficient_portfolio", env=env) as m:
            x = m.addMVar(shape=self.mu.shape, name="x")
            m.setObjective(
                self.mu @ x - 0.5 * gamma * (x @ (self.covariance @ x)), GRB.MAXIMIZE
            )

            investment = x.sum()

            b = m.addMVar(shape=self.mu.shape, vtype="B", name="trade_x")
            m.addConstr(x <= b)

            if max_trades is not None:
                m.addConstr(b.sum() <= max_trades)

            if fees_buy is not None:
                investment += b.sum() * fees_buy

            if min_buy_in is not None:
                m.addConstr(x >= min_buy_in * b, name="min_buy_in")

            m.addConstr(investment == 1, name="fully_invested")

            m.optimize()
            if m.Status == GRB.OPTIMAL:
                return self._convert_result(x.X)
