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
        self,
        gamma,
        max_trades=None,
        fees_buy=None,
        fees_sell=None,
        min_long=None,
        min_short=None,
        max_total_short=0.0,
    ):
        """
        Compute efficient portfolio for given paramters

        :param gamma: Risk aversion cofficient for balancing risk and return;
            the resulting objective functions is
            :math:`\mu^T x - 0.5 \gamma x^T \Sigma x`
        :type gamma: :class:`float` >= 0
        :param max_trades: Upper limit on the number of trades
        :type max_trades: :class:`int` >= 0
        :param fees_buy: Fixed-charge cost for each traded long position, relative
            to total portfolio value
        :type fees_buy: :class:`float` >= 0
        :param fees_sell: Fixed-charge cost for each traded short position, relative
            to total portfolio value
        :type fees_sell: :class:`float` >= 0
        :param min_long: Lower bound on the volume on a traded long position,
            relative to total portfolio value
        :type min_long: :class:`float` >= 0
        :param min_short: Lower bound on the volume on a traded short position,
            relative to total portfolio value
        :type min_long: :class:`float` >= 0
        :param max_total_short: Maximum total short positions, relative to
            total investment.
        :type max_total_short: :class:`float` >= 0
        """

        # max x' * mu + x' * Sigma * x
        # s.t.
        #      x = x_long - x_short  (x is split in positive/negative parts)
        #
        #      sum(x) + sum(b_long) * fees_buy + sum(b_short) * fees_sell = 1
        #                             (Fully invested, minus transaction fees)
        #
        #      x_long  <= M b_long    (force x_long to zero if not traded long)
        #                             (M >= 1 + max_total_short)
        #
        #      x_short <= M b_long    (force x_long to zero if not traded long)
        #                             (M >= max_total_short)
        #
        #      b_short + b_long <= 1  (Cannot go long and short at the same time)
        #
        #      sum(x_short) <= max_total_short   (Bound total leverage)
        #
        #      sum(b_short) + sum(b_long) <= max_trades  (Trade limit)
        #
        #      x_long >= min_long * b_long (minimum long position)
        #
        #      x_short >= min_short * b_short (minimum short position)
        #
        #    x free       (relative portfolio holdings)
        #    x_long  >= 0 (relative long holdings)
        #    x_short >= 0 (relative short holdings)
        #    b_long \in {0,1} (indicator variable for x_long)
        #    b_short \in {0,1} (indicator variable for x_short)

        with gp.Env() as env, gp.Model("efficient_portfolio", env=env) as m:
            # Portfolio vector x is split into long and short positions
            x = m.addMVar(shape=self.mu.shape, lb=-float("inf"), name="x")
            x_long = m.addMVar(shape=self.mu.shape, name="x_long")
            x_short = m.addMVar(shape=self.mu.shape, name="x_short")
            m.addConstr(x == x_long - x_short)

            # Binaries used to enforce VUB and minimum buy-in
            b_long = m.addMVar(shape=self.mu.shape, vtype="B", name="trade_long")
            b_short = m.addMVar(shape=self.mu.shape, vtype="B", name="trade_short")

            # Define VUB constraints for x_long and x_short.
            #
            # Going short by alpha means that each long position is upper
            # bounded by 1 + alpha, and each short position by alpha.
            # This is implied by the sum(x) == 1 constraint.
            m.addConstr(x_long <= (1.0 + max_total_short) * b_long)
            m.addConstr(x_short <= max_total_short * b_short)

            # A position can only by short or long, not both
            m.addConstr(b_long + b_short <= 1, name="long_or_short")

            # Bound total leverage
            m.addConstr(x_short.sum() <= max_total_short, name="total_short")

            investment = x.sum()

            if max_trades is not None:
                m.addConstr(b_long.sum() + b_short.sum() <= max_trades)

            if fees_buy is not None:
                investment += b_long.sum() * fees_buy

            if fees_sell is not None:
                investment += b_short.sum() * fees_sell

            if min_long is not None:
                m.addConstr(x_long >= min_long * b_long, name="min_long")

            if min_short is not None:
                m.addConstr(x_short >= min_short * b_short, name="min_short")

            m.addConstr(investment == 1, name="fully_invested")

            # Basic mean-variance weighted objective
            m.setObjective(
                self.mu @ x - 0.5 * gamma * (x @ (self.covariance @ x)), GRB.MAXIMIZE
            )

            m.optimize()
            if m.Status == GRB.OPTIMAL:
                return self._convert_result(x.X)
