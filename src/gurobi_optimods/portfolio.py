"""
Mean-Variance Portfolio
-----------------------
"""

import gurobipy as gp
from gurobipy import GRB
import numpy as np
import pandas as pd

from gurobi_optimods.utils import optimod


class MeanVariancePortfolio:
    """Optimal mean-variance portfolio solver.

    Instantiate an object of :class:`MeanVariancePortfolio` for given
    covariance matrix and return vector.  Use
    :meth:`MeanVariancePortfolio.efficient_portfolio` to solve for efficient
    portfolios with given parameters.

    Parameters
    ----------
    mu : 1-d :class:`np.ndarray`
        Vector of expected returns for each asset
    cov_matrix : 2-d :class:`np.ndarray`
        Covariance matrix :math:`\Sigma`
    cov_factors : tuple of :class:`np.ndarray`
        Covariance factors that constitute :math:`\Sigma`. Typically each
        element ``F`` of ``cov_matrix`` will either be a

            * n-by-k dense matrix, or a
            * n-by-n diagonal matrix.

        Each element ``F`` of ``cov_factors`` contributes the term ``F @ F.T``
        to :math:`\Sigma`; see also :ref:`factor models`

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
        rf_return=None,
        *,
        create_env,
    ):
        """Compute efficient portfolio for given parameters

        Parameters
        ----------

        gamma : :class:`float` >= 0
            Risk aversion cofficient for balancing risk and return; the
            resulting objective functions is
            :math:`\mu^T x - 0.5 \gamma x^T \Sigma x`
        max_trades : :class:`int` >= 0, optional
            Upper limit on the number of trades
        max_positions : :class:`int` >= 0, optional
            Upper limit on the number of open positions
        fees_buy : :class:`float` or :class:`np.ndarray` >= 0, optional
            Fixed-charge fee for each buy transaction, relative to total
            portfolio value
        fees_sell : :class:`float` or :class:`np.ndarray` >= 0, optional
            Fixed-charge fee for each sell transaction, relative to total
            portfolio value
        costs_buy : :class:`float` or :class:`np.ndarray` >= 0, optional
            Variable transaction costs for each buy transaction, relative to
            trade value
        costs_sell : :class:`float` or :class:`np.ndarray` >= 0, optional
            Variable transaction costs for each sell transaction, relative to
            trade value
        min_long : :class:`float` >= 0, optional
            Lower bound on the volume on a traded long position, relative to
            total portfolio value
        min_short : :class:`float` >= 0, optional
            Lower bound on the volume on a traded short position, relative to
            total portfolio value
        max_total_short : :class:`float` >= 0, optional
            Maximum total short positions, relative to total investment.
        initial_holdings : 1-d :class:`np.ndarray`, optional
            Initial portfolio holdings (sum needs to be <= 1)
        rf_return : :class:`float`, optional, default None
            Include a risk-free asset having return rate ``rf_return``.

        Returns
        -------
        mvp_result : dict
            A dict containing the efficient portfolio, along with auxiliary
            information:

            * ``mvp_result["x"]``: The portfolio vector :math:`x`
            * ``mvp_result["risk"]``: The estimated risk :math:`x^T \Sigma x`
              of the portfolio
            * ``mvp_result["return"]``: The estimated return :math:`\mu^T x` of
              the portfolio
            * ``mvp["x_rf"]`` relative investment in the risk-free asset.
              Present only if ``rf_return`` was non-None on input

            Some combinations of requested portfolio features may rule out
            **all** possible portfolios.  In this corner case the value
            ``None`` is returned.


        Notes
        -----
        Refer to :ref:`portfolio features` for a detailed discussion of all
        parameters.

        """

        fees_buy = self._homogenize_input(fees_buy)
        fees_sell = self._homogenize_input(fees_sell)
        costs_buy = self._homogenize_input(costs_buy)
        costs_sell = self._homogenize_input(costs_sell)
        initial_holdings = self._homogenize_input(initial_holdings)

        if initial_holdings is not None:
            if initial_holdings.sum() > 1.0:
                raise ValueError("Initial holding's sum must not exceed 1.0")
        else:
            initial_holdings = np.zeros(self.mu.shape)

        with create_env() as env, gp.Model("efficient_portfolio", env=env) as m:
            x, x_rf = self._populate_model(
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
                rf_return,
            )

            m.optimize()
            status = m.Status
            if status == GRB.OPTIMAL:
                x_vals = x.X
                x_rf_val = x_rf.X

        if status == GRB.OPTIMAL:
            return self._construct_result(x_vals, x_rf_val, rf_return)
        elif status in [GRB.INFEASIBLE, GRB.INF_OR_UNBD]:
            print("No portfolio satisfies the constraints!")
            return None
        else:
            return None

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
        rf_return,
    ):
        # max rf_return * x_rf + x' * mu - gamma * x' * cov_matrix * x
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
        #               + x_rf
        #      = 1
        #                             (fully invested, minus transaction costs and fees)
        #
        #      x_long  <= M b_long    (force x_long to zero if not traded long)
        #                             (M >= 1 + max_total_short)
        #
        #      b_short + b_long <= 1  (cannot go long and short at the same time)
        #      b_sell + b_buy <= 1  (cannot sell and buy at the same time)
        #
        #      sum(x_short) <= max_total_short   (bound total leverage)
        #
        #      sum(b_buy) + sum(b_sell) <= max_trades  (trade limit)
        #      sum(b_long) + sum(b_short) <= max_positions  (position limit)
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

        # Dummy variable for investment in risk-free asset,
        x_rf = m.addVar(lb=0.0, ub=0.0, name="x_rf")

        x_buy = m.addMVar(shape=self.mu.shape, name="x_buy")
        x_sell = m.addMVar(shape=self.mu.shape, name="x_sell")
        m.addConstr(x - initial_holdings == x_buy - x_sell)

        # Binaries used to enforce VUB and minimum position/trade size
        b_long = m.addMVar(shape=self.mu.shape, vtype="B", name="position_long")
        b_short = m.addMVar(shape=self.mu.shape, vtype="B", name="position_short")

        b_buy = m.addMVar(shape=self.mu.shape, vtype="B", name="trade_buy")
        b_sell = m.addMVar(shape=self.mu.shape, vtype="B", name="trade_sell")

        m.update()

        # Define VUB constraints for x_long and x_short.
        #
        # Going short by alpha means that each long position is upper
        # bounded by 1 + alpha, and each short position by alpha.
        # This is implied by the sum(x) == 1 constraint.
        m.addConstr(x_long <= (1.0 + max_total_short) * b_long)
        m.addConstr(x_short <= max_total_short * b_short)

        m.addConstr(x_buy <= (1.0 + max_total_short) * b_buy)
        m.addConstr(x_sell <= (1.0 + max_total_short) * b_sell)

        # A position/trade cannot be both short and long
        m.addConstr(b_long + b_short <= 1, name="long_or_short_position")
        m.addConstr(b_buy + b_sell <= 1, name="buy_or_sell")

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

        if rf_return is not None:
            x_rf.ub = 1.0
            investment += x_rf

        if min_long is not None:
            m.addConstr(x_buy >= min_long * b_buy, name="min_buy")

        if min_short is not None:
            m.addConstr(x_sell >= min_short * b_sell, name="min_sell")

        m.addConstr(investment == 1, name="fully_invested")

        if not isinstance(self.covariance, tuple):
            # Basic mean-variance weighted objective
            objexpr = self.mu @ x - 0.5 * gamma * x @ self.covariance @ x
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

        if rf_return is not None:
            objexpr += rf_return * x_rf

        m.setObjective(objexpr, GRB.MAXIMIZE)
        return (x, x_rf)

    def _construct_result(self, x, x_rf, rf_return):
        pf = dict()
        if self.resultType == "numpy":
            pf["x"] = x
        elif self.resultType == "pandas":
            pf["x"] = pd.Series(x, index=self.index)
        else:
            assert False

        pf["return"] = self.mu @ x

        if not isinstance(self.covariance, tuple):
            pf["risk"] = x @ self.covariance @ x
        else:
            pf["risk"] = 0.0
            for F in self.covariance:
                y = x @ F
                pf["risk"] += y @ y

        if rf_return is not None:
            pf["x_rf"] = x_rf
            pf["return"] += rf_return * x_rf

        return pf

    def _homogenize_input(self, input_data):
        # Check and unpack if input_data is a Series
        if isinstance(input_data, pd.Series):
            if self.index is not None:
                if any(self.index != input_data.index):
                    raise ValueError("Misaligned Series indexes: " + input_data.index)
            else:
                self.index = input_data.index
            input_data = input_data.to_numpy()

        return input_data
