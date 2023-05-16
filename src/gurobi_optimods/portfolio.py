import gurobipy as gp
from gurobipy import GRB
import numpy as np
import scipy.sparse as sp
import pandas as pd

from gurobi_optimods.utils import optimod


class MeanVariancePortfolio:
    """Optimal mean-variance portfolio solver.

    Instantiate an object of :class:`MeanVariancePortfolio` for given
    covariance matrix and return vector.  Use
    :meth:`MeanVariancePortfolio.efficient_portfolio` to solve for efficient
    portfolios with given parameters.

    :param mu: Return vector
    :type mu: 1-d :class:`np.ndarray`
    :param cov_matrix: Covariance matrix
    :type cov_matrix: 2-d :class:`np.ndarray`
    :param cov_factors: Covariance factors
    :type cov_factors: 2-d :class:`tuple` of :class:`np.ndarray`

    """

    def __init__(
        self,
        mu,
        cov_matrix=None,
        cov_factors=None,
    ):
        if cov_matrix is not None and cov_factors is not None:
            raise TypeError("Both cov_matrix and cov_factors given")

        if cov_matrix is not None:
            if isinstance(cov_matrix, pd.DataFrame):
                self.resultType = "pandas"
                self.index = cov_matrix.index
                self.covariance = cov_matrix.to_numpy()
            elif isinstance(cov_matrix, np.ndarray):
                self.covariance = cov_matrix
                self.resultType = "numpy"
            else:
                raise TypeError("Incompatible type of cov_matrix")
        elif cov_factors is not None:
            self.covariance = cov_factors
            self.index = None
        else:
            raise TypeError("No covariace data given")

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

    def _populate_model(
        self,
        m,
        gamma,
        max_trades,
        max_positions,
        fees_buy,
        fees_sell,
        costs_buy,
        costs_sell,
        min_long,
        min_short,
        max_total_short,
        initial_holdings,
    ):
        # max x' * mu + x' * cov_matrix * x
        # s.t.
        #      x = x_long - x_short  (x is split in positive/negative parts)
        #
        #      x_long = initial_holdings_long + x_long_buy - x_long_sell
        #      x_short == initial_holdings_short - x_short_buy + x_short_sell
        #      x - initial_holdings == x_buy - x_sell
        #
        #      sum(x)   + sum(b_buy) * fees_buy
        #               + sum(b_sell) * fees_sell
        #               + sum(x_buy) * costs_buy
        #               + sum(x_sell) * costs_sell
        #      = 1
        #                             (Fully invested, minus transaction costs and fees)
        #
        #      x_long  <= M b_long    (force x_long to zero if not traded long)
        #                             (M >= 1 + max_total_short)
        #
        #      b_short + b_long <= 1  (Cannot go long and short at the same time)
        #
        #      sum(x_short) <= max_total_short   (Bound total leverage)
        #
        #      sum(b_buy) + sum(b_sell) <= max_trades  (Trade limit)
        #      sum(b_long) + sum(b_short) <= max_positions
        #
        #      x_buy >= min_long * b_buy (minimum buy position)
        #
        #      x_sell >= min_short * b_sell (minimum sell position)
        #
        #    x free       (relative portfolio holdings)
        #    x_long  >= 0 (relative long holdings)
        #    x_short >= 0 (relative short holdings)
        #    x_buy >=0 (relative buy)
        #    x_sell >=0 (relative sell)
        #    b_long \in {0,1} (indicator variable for x_long)
        #    b_short \in {0,1} (indicator variable for x_short)
        #    b_buy \in {0,1} (indicator variable for x_buy)
        #    b_sell \in {0,1} (indicator variable for x_sell)

        # Portfolio vector x is split into long and short positions
        x = m.addMVar(shape=self.mu.shape, lb=-float("inf"), name="x")
        x_long = m.addMVar(shape=self.mu.shape, name="x_long")
        x_short = m.addMVar(shape=self.mu.shape, name="x_short")
        m.addConstr(x == x_long - x_short)

        x_buy = m.addMVar(shape=self.mu.shape, name="x_buy")
        x_sell = m.addMVar(shape=self.mu.shape, name="x_sell")
        m.addConstr(x - initial_holdings == x_buy - x_sell)

        # Binaries used to enforce VUB and minimum position/trade size
        b_long = m.addMVar(shape=self.mu.shape, vtype="B", name="position_long")
        b_short = m.addMVar(shape=self.mu.shape, vtype="B", name="position_short")

        b_buy = m.addMVar(shape=self.mu.shape, vtype="B", name="trade_buy")
        b_sell = m.addMVar(shape=self.mu.shape, vtype="B", name="trade_sell")

        # Define VUB constraints for x_long and x_short.
        #
        # Going short by alpha means that each long position is upper
        # bounded by 1 + alpha, and each short position by alpha.
        # This is implied by the sum(x) == 1 constraint.
        m.addConstr(x_long <= (1.0 + max_total_short) * b_long)
        m.addConstr(x_short <= max_total_short * b_short)

        m.addConstr(x_buy <= (1.0 + max_total_short) * b_buy)
        m.addConstr(x_sell <= (1.0 + max_total_short) * b_sell)

        # A position/trade can only by short or long, not both
        m.addConstr(b_long + b_short <= 1, name="long_or_short_position")
        m.addConstr(b_buy + b_sell <= 1, name="boy")

        # Bound total leverage
        m.addConstr(x_short.sum() <= max_total_short, name="total_short")

        investment = x.sum()

        if max_trades is not None:
            m.addConstr(b_buy.sum() + b_sell.sum() <= max_trades, name="max_trades")

        if max_positions is not None:
            m.addConstr(
                b_long.sum() + b_short.sum() <= max_positions, name="max_positions"
            )

        if fees_buy is not None:
            investment += (b_buy * fees_buy).sum()
        if fees_sell is not None:
            investment += (b_sell * fees_sell).sum()

        if costs_buy is not None:
            investment += (x_buy * costs_buy).sum()
        if costs_sell is not None:
            investment += (x_sell * costs_sell).sum()

        if min_long is not None:
            m.addConstr(x_buy >= min_long * b_buy, name="min_buy")

        if min_short is not None:
            m.addConstr(x_sell >= min_short * b_sell, name="min_sell")

        m.addConstr(investment == 1, name="fully_invested")

        if not isinstance(self.covariance, tuple):
            # Basic mean-variance weighted objective
            m.setObjective(
                self.mu @ x - 0.5 * gamma * x @ self.covariance @ x,
                GRB.MAXIMIZE,
            )
        else:
            # Idea:   We have given  Sigma =
            #
            #   factors[0] @ factors[0].T + ... + factors[l] @ factors[l].T
            #   =: F_0 @ F_0.T + ... + F_l @ F_l
            #
            # so that for each contributing term we can set
            #
            #   (x.T @ F_i) @ (F_i.T @ x) =: y_i @ y_i
            #
            # giving
            #
            #   min \sum_i gamma * y_i @ y_i
            #   s.t. y_i = F_i.T @ x  for all i

            objexpr = self.mu @ x

            for idx, F in enumerate(self.covariance):
                y = m.addMVar(F.shape[1], lb=-float("inf"), name=f"factor{idx:d}")
                m.addConstr(F.T @ x == y, name=f"link_factor{idx:d}_x")
                objexpr -= 0.5 * gamma * y @ y

            m.setObjective(objexpr, GRB.MAXIMIZE)

        return x

    @optimod()
    def efficient_portfolio(
        self,
        gamma,
        max_trades=None,
        max_positions=None,
        fees_buy=None,
        fees_sell=None,
        costs_buy=None,
        costs_sell=None,
        min_long=None,
        min_short=None,
        max_total_short=0.0,
        initial_holdings=None,
        gurobi_params=None,
        *,
        create_env,
    ):
        """
        Compute efficient portfolio for given paramters

        :param gamma: Risk aversion cofficient for balancing risk and return;
            the resulting objective functions is
            :math:`\mu^T x - 0.5 \gamma x^T \cov_matrix x`
        :type gamma: :class:`float` >= 0
        :param max_trades: Upper limit on the number of trades
        :type max_trades: :class:`int` >= 0
        :param max_positions: Upper limit on the number of open positions
        :type max_positions: :class:`int` >= 0
        :param fees_buy: Fixed-charge fee for each buy transaction, relative
            to total portfolio value
        :type fees_buy: :class:`float` or :class:`np.ndarray` >= 0
        :param fees_sell: Fixed-charge fee for each sell transaction, relative
            to total portfolio value
        :type fees_buy: :class:`float` or :class:`np.ndarray` >= 0
        :param costs_buy: Variable transaction costs for each buy transaction, relative
            to trade value
        :type fees_buy: :class:`float` or :class:`np.ndarray` >= 0
        :param costs_sell: Variable transaction costs for each sell transaction, relative
            to trade value
        :type fees_buy: :class:`float` or :class:`np.ndarray` >= 0
        :param min_long: Lower bound on the volume on a traded long position,
            relative to total portfolio value
        :type min_long: :class:`float` >= 0
        :param min_short: Lower bound on the volume on a traded short position,
            relative to total portfolio value
        :type min_long: :class:`float` >= 0
        :param max_total_short: Maximum total short positions, relative to
            total investment.
        :type max_total_short: :class:`float` >= 0
        :param max_total_short: Maximum total short positions, relative to
            total investment.
        :type max_total_short: :class:`float` >= 0
        :param initial_holdings: Initial portfolio holdings (sum needs to be <= 1)
        :type initial_holdings: 1-d :class:`np.ndarray`
        :param gurobi_params: Gurobi parameters to be passed to the solver
        :type gurobi_params: class:`dict`
        :param silent: silent=True suppresses all console output (defaults to False)
        :type silent: bool

        Refer to the Section :ref:`portfolio features` for a detailed discussion
        of these parameters.
        """

        if isinstance(initial_holdings, pd.Series):
            initial_holdings = initial_holdings.to_numpy()

        if isinstance(fees_buy, pd.Series):
            fees_buy = fees_buy.to_numpy()
        if isinstance(fees_sell, pd.Series):
            fees_sell = fees_sell.to_numpy()
        if isinstance(costs_buy, pd.Series):
            costs_buy = costs_buy.to_numpy()
        if isinstance(costs_sell, pd.Series):
            costs_sell = costs_sell.to_numpy()

        if initial_holdings is not None:
            if initial_holdings.sum() > 1.0:
                raise ValueError("Initial holding's sum must not exceed 1.0")
        else:
            initial_holdings = np.zeros(self.mu.shape)

        with create_env(params=gurobi_params) as env, gp.Model(
            "efficient_portfolio", env=env
        ) as m:
            x = self._populate_model(
                m,
                gamma,
                max_trades,
                max_positions,
                fees_buy,
                fees_sell,
                costs_buy,
                costs_sell,
                min_long,
                min_short,
                max_total_short,
                initial_holdings,
            )

            m.optimize()
            status = m.Status
            if status == GRB.OPTIMAL:
                xvals = x.X

        if status == GRB.OPTIMAL:
            return self._convert_result(xvals)
        else:
            return None
