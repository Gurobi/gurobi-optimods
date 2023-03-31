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


# very simple Markowitz type mean-variance portfolio
class MeanVariancePortfolio:
    def __init__(
        self,
        Sigma,
        mu,
    ) -> None:
        if isinstance(Sigma, pd.DataFrame):
            self.covariance = Sigma.to_numpy()
        elif isinstance(Sigma, np.ndarray):
            self.covariance = Sigma
        else:
            raise TypeError("Incompatible type of Sigma")

        if isinstance(mu, pd.Series):
            self.mu = mu.to_numpy()
        elif isinstance(mu, np.ndarray):
            self.mu = mu
        else:
            raise TypeError("Incompatible type of mu")

    def minimize_risk(self, expected_return):
        try:
            with gp.Env() as env, gp.Model("min_risk", env=env) as m:
                x = m.addMVar(shape=self.mu.shape, name="x")

                m.addConstr(x.sum() == 1, name="fully_invested")
                m.addConstr(self.mu @ x >= expected_return, name="expected_return")
                m.setObjective(x @ self.covariance @ x)

                m.optimize()

                if m.Status == GRB.OPTIMAL:
                    return x.X

        except gp.GurobiError as e:
            print("Error code " + str(e.errno) + ": " + str(e))
        except AttributeError:
            print("Encountered an attribute error")

    def maximize_return(self, max_risk):
        try:
            with gp.Env() as env, gp.Model("max_return", env=env) as m:
                x = m.addMVar(shape=self.mu.shape, name="x")
                m.addConstr(x.sum() == 1, name="fully_invested")
                m.addConstr(x @ self.covariance @ x <= max_risk, name="max_risk")
                m.setObjective(self.mu @ x, GRB.MAXIMIZE)

                m.optimize()

                if m.Status == GRB.OPTIMAL:
                    return x.X

        except gp.GurobiError as e:
            print("Error code " + str(e.errno) + ": " + str(e))
        except AttributeError:
            print("Encountered an attribute error")

    # trade-off between minimizing risk and maximizing return for a given risk-aversion coefficient gamma
    def efficient_portfolio(self, gamma):
        try:
            with gp.Env() as env, gp.Model("efficient_portfolio", env=env) as m:
                x = m.addMVar(shape=self.mu.shape, name="x")
                m.addConstr(x.sum() == 1, name="fully_invested")
                m.setObjective(
                    self.mu @ x - 0.5 * gamma * x @ self.covariance @ x, GRB.MAXIMIZE
                )

                m.optimize()
                if m.Status == GRB.OPTIMAL:
                    return x.X

        except gp.GurobiError as e:
            print("Error code " + str(e.errno) + ": " + str(e))
        except AttributeError:
            print("Encountered an attribute error")
