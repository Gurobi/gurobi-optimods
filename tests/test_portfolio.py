from contextlib import redirect_stdout
import io
import unittest
import numpy as np
import pandas as pd
from numpy.testing import assert_allclose, assert_array_equal

from gurobi_optimods.datasets import load_portfolio
from gurobi_optimods.portfolio import MeanVariancePortfolio, PortfolioResult


class TestMVPBasic(unittest.TestCase):
    # All tests that focus on construction, data checks etc. go here

    def test_datasets(self):
        data = load_portfolio()
        self.assertEqual(
            set(data.keys()),
            {"AA", "BB", "CC", "DD", "EE", "FF", "GG", "HH", "II", "JJ"},
        )

    def test_example_data(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()

        mvp = MeanVariancePortfolio(mu, cov_matrix)
        pf = mvp.efficient_portfolio(0.5)

    def test_init_0(self):
        # Specifying neither cov_matrix nor cov_factors is disallowed
        data = load_portfolio()
        mu = data.mean()

        with self.assertRaises(TypeError):
            mvp = MeanVariancePortfolio(mu)

        with self.assertRaises(TypeError):
            mvp = MeanVariancePortfolio(mu, cov_matrix=None)

        with self.assertRaises(TypeError):
            mvp = MeanVariancePortfolio(mu, cov_factors=None)

        with self.assertRaises(TypeError):
            mvp = MeanVariancePortfolio(mu, cov_matrix=None, cov_factors=None)

    def test_init_1(self):
        # Specifying both cov_matrix and cov_factors is disallowed
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        L = np.linalg.cholesky(cov_matrix)

        with self.assertRaises(TypeError):
            mvp = MeanVariancePortfolio(mu, cov_matrix=cov_matrix, cov_factors=(L,))

    def test_outputflag(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()

        mvp = MeanVariancePortfolio(mu, cov_matrix)

        with redirect_stdout(io.StringIO()) as console:
            pf = mvp.efficient_portfolio(0.5, solver_params={})
        consoleContent = console.getvalue()
        self.assertIn("Gurobi Optimizer", consoleContent)

        solver_params = {"OutputFlag": 0}
        with redirect_stdout(io.StringIO()) as console:
            pf = mvp.efficient_portfolio(0.5, solver_params=solver_params)
        consoleContent = console.getvalue()
        self.assertEqual(consoleContent, "")

    def test_non_verbose(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()

        mvp = MeanVariancePortfolio(mu, cov_matrix)

        with redirect_stdout(io.StringIO()) as console:
            pf = mvp.efficient_portfolio(0.5, solver_params={})
        consoleContent = console.getvalue()
        self.assertIn("Gurobi Optimizer", consoleContent)

        with redirect_stdout(io.StringIO()) as console:
            pf = mvp.efficient_portfolio(0.5, verbose=False)
        consoleContent = console.getvalue()
        self.assertEqual(consoleContent, "")

    def test_logfile(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()

        mvp = MeanVariancePortfolio(mu, cov_matrix)

        pf = mvp.efficient_portfolio(0.5, logfile="more_bananas.log")
        with open("more_bananas.log", "rt") as fh:
            content = fh.read()
        self.assertIn("Gurobi Optimizer", content)

    def test_params(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()

        mvp = MeanVariancePortfolio(mu, cov_matrix)

        solver_params = {"MIPFocus": 1}
        with redirect_stdout(io.StringIO()) as console:
            pf = mvp.efficient_portfolio(0.5, solver_params=solver_params)
        consoleContent = console.getvalue()
        self.assertIn("Set parameter MIPFocus to value 1", consoleContent)


class TestMVPFeatures(unittest.TestCase):
    # All tests that focus portfolio feature correctness go here

    def test_two_assets_x(self):
        cov_matrix = np.array([[3, 0.5], [0.5, 2]])
        mu = np.array([1, -0.1])
        mvp = MeanVariancePortfolio(mu, cov_matrix)
        pf = mvp.efficient_portfolio(0.5)
        self.assertTrue(hasattr(pf, "x"))
        assert_allclose(pf.x, [0.925, 0.075], atol=1e-6)

    def test_two_assets_return(self):
        cov_matrix = np.array([[3, 0.5], [0.5, 2]])
        mu = np.array([1, -0.1])
        mvp = MeanVariancePortfolio(mu, cov_matrix)
        pf = mvp.efficient_portfolio(0.5)
        self.assertTrue(hasattr(pf, "ret"))
        self.assertAlmostEqual(pf.ret, pf.x @ mu)

    def test_two_assets_risk(self):
        cov_matrix = np.array([[3, 0.5], [0.5, 2]])
        mu = np.array([1, -0.1])
        mvp = MeanVariancePortfolio(mu, cov_matrix)
        pf = mvp.efficient_portfolio(0.5)
        self.assertTrue(hasattr(pf, "risk"))
        self.assertAlmostEqual(pf.risk, pf.x @ cov_matrix @ pf.x)

    def test_example_data_result(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0

        mvp = MeanVariancePortfolio(mu, cov_matrix)
        pf = mvp.efficient_portfolio(gamma)
        self.assertIsInstance(pf, PortfolioResult)
        self.assertTrue(hasattr(pf, "x"))
        self.assertTrue(hasattr(pf, "ret"))
        self.assertTrue(hasattr(pf, "risk"))
        self.assertTrue(hasattr(pf, "x_rf"))

        self.assertAlmostEqual(pf.ret, pf.x @ mu)
        self.assertAlmostEqual(pf.risk, pf.x @ cov_matrix @ pf.x)
        self.assertIsNone(pf.x_rf)

    def test_number_of_trades(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0

        mvp = MeanVariancePortfolio(mu, cov_matrix)
        x_unconstrained = mvp.efficient_portfolio(gamma).x
        self.assertGreater((x_unconstrained > 1e-4).sum(), 3)

        x_3 = mvp.efficient_portfolio(gamma, max_trades=3).x
        self.assertLessEqual((x_3 > 1e-4).sum(), 3)

    def test_transaction_fees_long(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0
        fees_buy = 1e-4

        mvp = MeanVariancePortfolio(mu, cov_matrix)
        x = mvp.efficient_portfolio(gamma, fees_buy=fees_buy).x
        n_trades = (x > 1e-4).sum()
        self.assertLessEqual(x.sum(), 1 - n_trades * fees_buy)

    def test_transaction_fees_short(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0
        fees_sell = 1e-4
        mvp = MeanVariancePortfolio(mu, cov_matrix)

        # Ensure that we go short somewhere
        x = mvp.efficient_portfolio(gamma, max_total_short=0.3).x
        n_trades = (x < 0).sum()
        self.assertGreater(n_trades, 0)
        self.assertAlmostEqual(x.sum(), 1)

        # Ensure that transaction fees are paid out of the portfolio
        x = mvp.efficient_portfolio(gamma, max_total_short=0.3, fees_sell=fees_sell).x
        n_trades = (x < 0).sum()
        self.assertGreater(n_trades, 0)
        self.assertLessEqual(x.sum(), 1 - n_trades * fees_sell + 1e-8)

    def test_min_long(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0

        mvp = MeanVariancePortfolio(mu, cov_matrix)

        # Ensure that we _have_ a small trade w/o minimum buy in
        x = mvp.efficient_portfolio(gamma, min_long=0.0).x
        small_trades = (x > 1e-6) & (x < (0.03 - 1e-6))
        self.assertGreater(small_trades.sum(), 0)

        # Ensure that we _don't_ have a small trade w minimum buy in
        x = mvp.efficient_portfolio(gamma, min_long=0.03).x
        small_trades = (x > 1e-6) & (x < (0.03 - 1e-6))
        self.assertEqual(small_trades.sum(), 0)

    def test_max_total_short(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0

        mvp = MeanVariancePortfolio(mu, cov_matrix)

        # Ensure that by default we don't go short
        x = mvp.efficient_portfolio(gamma).x
        self.assertEqual((x < 0).sum(), 0)

        # Ensure that we take advantage of leverage
        x = mvp.efficient_portfolio(gamma, max_total_short=0.1).x
        self.assertGreaterEqual(x[x < 0].sum(), -0.1 - 1e-6)
        self.assertLess(x[x < 0].sum(), -1e-3)
        self.assertAlmostEqual(x.sum(), 1.0)
        self.assertAlmostEqual(np.abs(x).sum(), 1.0 + 2 * 0.1)

    def test_min_short_0(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0

        mvp = MeanVariancePortfolio(mu, cov_matrix)

        # Ensure that by default we don't go short
        x = mvp.efficient_portfolio(gamma).x
        self.assertEqual((x < 0).sum(), 0)

        # Adding a min_short constraint doesn't change a thing
        x_other = mvp.efficient_portfolio(gamma, min_short=0.01).x
        assert_allclose(x.to_numpy(), x_other.to_numpy())

    def test_min_short_1(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0

        mvp = MeanVariancePortfolio(mu, cov_matrix)
        x = mvp.efficient_portfolio(gamma, max_total_short=0.5).x

        # Ensure that we do have tiny short positions
        small_trades = (x < -1e-6) & (x > (-0.05 + 1e-6))
        self.assertGreater(small_trades.sum(), 0)

        # Ensure that we still have short positions but beyond the threshold
        x = mvp.efficient_portfolio(gamma, max_total_short=0.5, min_short=0.05).x
        self.assertGreater((x <= -0.05).sum(), 0)
        small_trades = (x < -1e-6) & (x > (-0.05 + 1e-6))
        self.assertEqual(small_trades.sum(), 0)

    def test_start_portfolio_empty(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0

        mvp = MeanVariancePortfolio(mu, cov_matrix)
        x0 = np.zeros(mu.shape)
        x_with = mvp.efficient_portfolio(gamma, initial_holdings=x0).x
        x_without = mvp.efficient_portfolio(gamma, initial_holdings=None).x
        assert_array_equal(x_with, x_without)

    def test_start_portfolio_invalid(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0

        mvp = MeanVariancePortfolio(mu, cov_matrix)
        x0 = np.ones(mu.shape)

        with self.assertRaises(ValueError):
            x_with = mvp.efficient_portfolio(gamma, initial_holdings=x0)

    def test_start_portfolio_without_restrictions(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0

        # If there are no additional restrictions, the resulting portfolio
        # should be the same as without initial holdings sunk-cost-fallacy
        mvp = MeanVariancePortfolio(mu, cov_matrix)
        x0 = 1.0 / mu.size * np.ones(mu.size)
        x_with = mvp.efficient_portfolio(gamma, initial_holdings=x0).x
        x_without = mvp.efficient_portfolio(gamma, initial_holdings=None).x
        assert_allclose(x_with, x_without, atol=1e-6)

    def test_start_portfolio_not_fully_invested(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0

        # If there are no additional restrictions, the resulting portfolio
        # should be the same as without initial holdings sunk-cost-fallacy
        mvp = MeanVariancePortfolio(mu, cov_matrix)
        x0 = 0.5 / mu.size * np.ones(mu.size)
        x_with = mvp.efficient_portfolio(gamma, initial_holdings=x0).x
        x_without = mvp.efficient_portfolio(gamma, initial_holdings=None).x
        assert_allclose(x_with, x_without, atol=1e-6)

    def test_start_portfolio_limit_trades(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0

        mvp = MeanVariancePortfolio(mu, cov_matrix)
        x0 = 1.0 / mu.size * np.ones(mu.size)
        # Ensure that we use more than 3 trades if there is no limit
        x1_nomaxtrades = mvp.efficient_portfolio(gamma, initial_holdings=x0).x
        trades = ((x1_nomaxtrades - x0) > 1e-4) | ((x1_nomaxtrades - x0) < -1e-4)
        self.assertGreater(trades.sum(), 3)

        # Ensure that we only use 3 trades with the limit
        x1_maxtrades = mvp.efficient_portfolio(
            gamma, initial_holdings=x0, max_trades=3
        ).x
        trades = ((x1_maxtrades - x0) > 1e-4) | ((x1_maxtrades - x0) < -1e-4)
        self.assertLessEqual(trades.sum(), 3)

    def test_start_portfolio_max_total_short_max_trades(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0

        mvp = MeanVariancePortfolio(mu, cov_matrix)
        x0 = 1.0 / mu.size * np.ones(mu.size)

        # Ensure that we take advantage of leverage
        x = mvp.efficient_portfolio(
            gamma, initial_holdings=x0, max_trades=6, max_total_short=0.1
        ).x

        self.assertGreaterEqual((x[x < 0]).sum(), -0.1 - 1e-6)
        self.assertLess(x[x < 0].sum(), -1e-3)
        self.assertAlmostEqual(x.sum(), 1.0)

        trades = ((x - x0) > 1e-4) | ((x - x0) < -1e-4)
        self.assertLessEqual(trades.sum(), 6)

    def test_start_portfolio_min_buy(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0

        mvp = MeanVariancePortfolio(mu, cov_matrix)
        x0 = np.array([0.4, 0, 0, 0.2, 0, 0.1, 0.03, 0.2, 0.07, 0])
        x0 /= x0.sum()  # To avoid 1+eps results due to rounding
        # Ensure that we _do_ have a small trade w/o minimum buy
        x = mvp.efficient_portfolio(gamma, min_long=0.0, initial_holdings=x0).x
        trades = x - x0
        small_trades = (trades > 1e-6) & (trades < (0.03 - 1e-6))
        self.assertGreater(small_trades.sum(), 0)

        # Ensure that we _don't_ have a small trade w/ minimum buy
        x = mvp.efficient_portfolio(gamma, min_long=0.03, initial_holdings=x0).x
        trades = x - x0
        small_trades = (trades > 1e-6) & (trades < (0.03 - 1e-6))
        self.assertEqual(small_trades.sum(), 0)

    def test_start_portfolio_no_fees(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0
        fees = 1e-4

        mvp = MeanVariancePortfolio(mu, cov_matrix)

        # Determine efficient portfolio without fees
        x0 = mvp.efficient_portfolio(gamma, max_total_short=0.1).x

        # Ensure that this is not changed when we compute the portfolio again
        # with x0 as start and fees.  We need min_long and min_short to avoid
        # that paying fees is seen as a "risk-free loss" which might improve
        # our objective
        x = mvp.efficient_portfolio(
            gamma,
            max_total_short=0.1,
            initial_holdings=x0.to_numpy(),
            fees_buy=fees,
            fees_sell=fees,
            min_long=0.02,
            min_short=0.02,
        ).x
        assert_allclose(x0, x, atol=1e-6)

    def test_start_portfolio_fees_buy(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0
        fees = 1e-4

        mvp = MeanVariancePortfolio(mu, cov_matrix)
        x0 = 1.0 / mu.size * np.ones(mu.size)

        x = mvp.efficient_portfolio(
            gamma,
            max_total_short=0.1,
            initial_holdings=x0,
            fees_buy=fees,
            min_long=0.02,
            min_short=0.02,
        ).x
        buy_trades = (x - x0) > 1e-4
        # Ensure that there have been buy trades
        self.assertGreater(buy_trades.sum(), 0)
        self.assertLessEqual(x.sum(), 1 - buy_trades.sum() * fees + 1e-6)

    def test_start_portfolio_fees_sell(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0
        fees = 1e-4

        mvp = MeanVariancePortfolio(mu, cov_matrix)
        x0 = 1.0 / mu.size * np.ones(mu.size)

        x = mvp.efficient_portfolio(
            gamma,
            max_total_short=0.1,
            initial_holdings=x0,
            fees_sell=fees,
            min_long=0.02,
            min_short=0.02,
        ).x
        sell_trades = (x0 - x) > 1e-4
        # Ensure that there have been sell trades
        self.assertGreater(sell_trades.sum(), 0)
        self.assertLessEqual(x.sum(), 1 - sell_trades.sum() * fees + 1e-6)

    def test_max_positions(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0

        mvp = MeanVariancePortfolio(mu, cov_matrix)
        x = mvp.efficient_portfolio(gamma, max_positions=3).x
        self.assertLessEqual((x > 1e-4).sum(), 3)

    def test_max_positions_start_portfolio(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0

        mvp = MeanVariancePortfolio(mu, cov_matrix)
        x0 = 1.0 / mu.size * np.ones(mu.size)

        x = mvp.efficient_portfolio(
            gamma,
            initial_holdings=x0,
            max_positions=6,
            max_trades=5,
        ).x

        self.assertLessEqual((x > 1e-4).sum(), 6)
        # In order to satisfy max_positions=6 with max_trades=5,
        # we need to sell 4 positions and invest the surplus into one of the remaining open positions.
        self.assertEqual((x >= (5.0 / mu.size - 1e-6)).sum(), 1)

    def test_max_positions_infeasible(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0

        mvp = MeanVariancePortfolio(mu, cov_matrix)
        x0 = 1.0 / mu.size * np.ones(mu.size)

        pf = mvp.efficient_portfolio(
            gamma,
            initial_holdings=x0,
            max_positions=6,
            max_trades=3,
        )
        # This needs to be infeasible.
        self.assertIsNone(pf)

    def test_transaction_costs_long(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0
        costs = 0.0025

        mvp = MeanVariancePortfolio(mu, cov_matrix)
        x = mvp.efficient_portfolio(gamma, costs_buy=costs).x
        value_trades = x.sum()
        self.assertAlmostEqual(x.sum(), 1 - costs * value_trades, delta=1e-6)

    def test_transaction_fees_costs_long(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0
        costs = 0.0025
        fees = 1e-4

        mvp = MeanVariancePortfolio(mu, cov_matrix)
        x = mvp.efficient_portfolio(gamma, costs_buy=costs, fees_buy=fees).x
        n_trades = (x > 1e-4).sum()
        value_trades = x.sum()
        self.assertAlmostEqual(
            x.sum(), 1 - n_trades * fees - value_trades * costs, delta=1e-3
        )

    def test_transaction_fees_costs_short(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0
        fees = 1e-4
        costs = 0.0025
        mvp = MeanVariancePortfolio(mu, cov_matrix)

        x = mvp.efficient_portfolio(
            gamma,
            max_total_short=0.3,
            fees_sell=fees,
            costs_sell=costs,
        ).x
        n_trades = (x < 0).sum()
        value_trades = -x[x < 0].sum()

        # Ensure that we go short somewhere
        self.assertGreater(n_trades, 0)
        # Ensure that transaction fees and costs are paid out of the portfolio
        self.assertAlmostEqual(x.sum(), 1 - n_trades * fees - value_trades * costs)

    def test_transaction_fees_costs_all(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0
        fees = 1e-4
        costs = 0.0025
        mvp = MeanVariancePortfolio(mu, cov_matrix)

        x = mvp.efficient_portfolio(
            gamma,
            max_total_short=0.3,
            fees_buy=fees,
            fees_sell=fees,
            costs_buy=costs,
            costs_sell=costs,
        ).x
        trades = (x < -1e-6) | (x > 1e-6)
        n_trades = trades.sum()
        long_trades_value = x[x > 0].sum()
        short_trades_value = -x[x < 0].sum()
        self.assertAlmostEqual(
            x.sum(),
            1 - n_trades * fees - (long_trades_value + short_trades_value) * costs,
        )

    def test_start_portfolio_transaction_fees_costs_all(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0
        fees = 1e-4
        costs = 0.0025
        mvp = MeanVariancePortfolio(mu, cov_matrix)
        x0 = 1.0 / mu.size * np.ones(mu.size)

        x = mvp.efficient_portfolio(
            gamma,
            max_total_short=0.1,
            initial_holdings=x0,
            fees_buy=fees,
            fees_sell=fees,
            costs_buy=costs,
            costs_sell=costs,
            min_long=0.02,
            min_short=0.02,
            max_trades=6,
        ).x
        trades = x - x0
        buy_trades = trades > 1e-4
        sell_trades = trades < -1e-4
        n_trades = buy_trades.sum() + sell_trades.sum()
        buy_trades_value = trades[buy_trades].sum()
        sell_trades_value = -trades[sell_trades].sum()
        trades_value = buy_trades_value + sell_trades_value
        self.assertAlmostEqual(x.sum(), 1 - n_trades * fees - trades_value * costs)

    def test_start_portfolio_max_total_short(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0
        mvp = MeanVariancePortfolio(mu, cov_matrix)

        # Starting portfolio that is short in one position, matching the
        # allowed leverage in the optimization
        x0 = np.zeros(mu.size)
        x0[2] = -0.1

        x = mvp.efficient_portfolio(
            gamma,
            initial_holdings=x0,
            min_short=0.01,
            min_long=0.015,
            max_total_short=0.1,
            fees_buy=0.001,
            fees_sell=0.002,
        ).x

        self.assertAlmostEqual(x[x < 0].sum(), -0.1)

    def test_risk_factors_equivalent_0(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        L = np.linalg.cholesky(cov_matrix)
        mu = data.mean()
        gamma = 100.0

        # Two equivalent setups: Sigma itself, and its Cholesky factor
        mvp_factors = MeanVariancePortfolio(mu, None, cov_factors=(L,))
        mvp_Sigma = MeanVariancePortfolio(mu, cov_matrix)

        x_factors = mvp_factors.efficient_portfolio(gamma).x
        x_Sigma = mvp_Sigma.efficient_portfolio(gamma).x

        assert_allclose(x_factors, x_Sigma, atol=1e-6)

    def test_risk_factors_auxdata_0(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        L = np.linalg.cholesky(cov_matrix)
        mu = data.mean()
        gamma = 100.0

        mvp = MeanVariancePortfolio(mu, None, cov_factors=(L,))
        pf = mvp.efficient_portfolio(gamma)

        self.assertAlmostEqual(pf.ret, mu.to_numpy() @ pf.x)
        self.assertAlmostEqual(pf.risk, cov_matrix.to_numpy() @ pf.x @ pf.x)

    def test_risk_factors_equivalent_1(self):
        # Use case: Data that would emerge from a Single factor model
        # R = alpha + beta * R_m + e.
        #
        # Letting s = beta * R_m, the resulting covariance matrix would then
        # become Sigma = s @ s.T + diag(...) where the second term captures the
        # variances of the specific returns.
        #
        # We don't bother simulating any such data, we just make up the data.

        n = 3

        mu = 0.2 + (0.2 * np.random.rand(3) - 0.1)
        s = (0.1 + 0.2 * np.random.rand(3)).reshape((3, 1))
        d = 0.05 + 0.3 * np.random.rand(3)
        D = np.diag(d)
        cov_matrix = s @ s.T + np.diag(d**2)

        # Two equivalent setups: Sigma itself, and its Cholesky factor
        gamma = 20
        mvp_factors = MeanVariancePortfolio(mu, None, cov_factors=(s, D))
        mvp_Sigma = MeanVariancePortfolio(mu, cov_matrix)

        x_factors = mvp_factors.efficient_portfolio(gamma).x
        x_Sigma = mvp_Sigma.efficient_portfolio(gamma).x

        assert_allclose(x_factors, x_Sigma, atol=1e-5)

    def test_risk_factors_auxdata_1(self):
        n = 3
        mu = 0.2 + (0.2 * np.random.rand(3) - 0.1)
        s = (0.1 + 0.2 * np.random.rand(3)).reshape((3, 1))
        d = 0.05 + 0.3 * np.random.rand(3)
        D = np.diag(d)
        cov_matrix = s @ s.T + np.diag(d**2)

        gamma = 20
        mvp = MeanVariancePortfolio(mu, None, cov_factors=(s, D))
        pf = mvp.efficient_portfolio(gamma)

        self.assertAlmostEqual(pf.ret, mu @ pf.x)
        self.assertAlmostEqual(pf.risk, cov_matrix @ pf.x @ pf.x)

    def test_costs_per_asset_long(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0
        mvp = MeanVariancePortfolio(mu, cov_matrix)

        costs = np.array(
            [0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.0025, 0.0025, 0.0025, 0.0025]
        )
        x = mvp.efficient_portfolio(gamma, costs_buy=costs).x

        self.assertAlmostEqual(x.sum(), 1 - (costs * x.to_numpy()).sum(), delta=1e-6)

    def test_costs_per_asset_long_random(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0
        mvp = MeanVariancePortfolio(mu, cov_matrix)

        costs = np.random.rand(mu.size) * 0.05
        x = mvp.efficient_portfolio(gamma, costs_buy=costs).x

        self.assertAlmostEqual(x.sum(), 1 - (costs * x.to_numpy()).sum(), delta=1e-6)

    def test_costs_per_asset_long_and_short(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0
        mvp = MeanVariancePortfolio(mu, cov_matrix)

        costs_buy = np.array(
            [0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.0025, 0.0025, 0.0025, 0.0025]
        )
        costs_sell = np.array(
            [0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.0025, 0.0025, 0.0025, 0.0025]
        )
        x = mvp.efficient_portfolio(
            gamma, costs_buy=costs_buy, costs_sell=costs_sell, max_total_short=0.1
        ).x

        df_buy = pd.Series(costs_buy, index=mu.index)
        df_sell = pd.Series(costs_sell, index=mu.index)
        total_costs_buy = (df_buy[x > 1e-6] * x[x > 1e-6]).sum()
        total_costs_sell = -(df_sell[x < -1e-6] * x[x < -1e-6]).sum()
        self.assertAlmostEqual(
            x.sum(), 1 - total_costs_buy - total_costs_sell, delta=1e-6
        )

    def test_costs_per_asset_long_and_short_random(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0
        mvp = MeanVariancePortfolio(mu, cov_matrix)

        # costs = np.array([0.00832292, 0.00295583, 0.00353887, 0.0044858 , 0.00337871, 0.00924592, 0.00339121, 0.00391734, 0.00773179, 0.00425846])
        costs_buy = np.random.rand(mu.size) * 0.05
        costs_sell = np.random.rand(mu.size) * 0.05
        x = mvp.efficient_portfolio(
            gamma, costs_buy=costs_buy, costs_sell=costs_sell, max_total_short=0.1
        ).x

        df_buy = pd.Series(costs_buy, index=mu.index)
        df_sell = pd.Series(costs_sell, index=mu.index)
        total_costs_buy = (df_buy[x > 1e-6] * x[x > 1e-6]).sum()
        total_costs_sell = -(df_sell[x < -1e-6] * x[x < -1e-6]).sum()
        self.assertAlmostEqual(
            x.sum(), 1 - total_costs_buy - total_costs_sell, delta=1e-6
        )

    def test_fees_per_asset_long(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0
        mvp = MeanVariancePortfolio(mu, cov_matrix)

        fees = np.array(
            [0.01, 0.003, 0.01, 0.004, 0.001, 0.0002, 0.002, 0.006, 0.005, 0.004]
        )
        x = mvp.efficient_portfolio(gamma, fees_buy=fees, max_trades=4).x
        df_buy = pd.Series(fees, index=mu.index)
        total_fees = df_buy[x > 1e-6].sum()

        self.assertAlmostEqual(x.sum(), 1 - total_fees, delta=1e-6)

    def test_fees_per_asset_long_random(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0
        mvp = MeanVariancePortfolio(mu, cov_matrix)

        fees = np.random.rand(10) * 0.01
        x = mvp.efficient_portfolio(gamma, fees_buy=fees, max_trades=4).x
        df_buy = pd.Series(fees, index=mu.index)
        total_fees = df_buy[x > 1e-6].sum()

        self.assertAlmostEqual(x.sum(), 1 - total_fees, delta=1e-6)

    def test_fees_per_asset_long_and_short_random(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0
        mvp = MeanVariancePortfolio(mu, cov_matrix)

        fees_buy = np.random.rand(10) * 0.01
        fees_sell = np.random.rand(10) * 0.01
        x = mvp.efficient_portfolio(
            gamma,
            fees_buy=fees_buy,
            fees_sell=fees_sell,
            max_trades=4,
            max_total_short=0.1,
        ).x
        df_buy = pd.Series(fees_buy, index=mu.index)
        df_sell = pd.Series(fees_sell, index=mu.index)
        total_fees = df_buy[x > 1e-6].sum() + df_sell[x < -1e-6].sum()

        self.assertAlmostEqual(x.sum(), 1 - total_fees, delta=1e-6)

    def test_input_wrong_index(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0
        mvp = MeanVariancePortfolio(mu, cov_matrix)

        fees_buy = pd.Series(index=mu.index[::-1], data=np.linspace(0.0, 0.1, mu.size))
        with self.assertRaises(ValueError):
            _ = mvp.efficient_portfolio(gamma, fees_buy=fees_buy)

        fees_sell = pd.Series(index=mu.index[::-1], data=np.linspace(0.0, 0.1, mu.size))
        with self.assertRaises(ValueError):
            _ = mvp.efficient_portfolio(gamma, fees_sell=fees_sell)

        costs_buy = pd.Series(index=mu.index[::-1], data=np.linspace(0.0, 0.1, mu.size))
        with self.assertRaises(ValueError):
            _ = mvp.efficient_portfolio(gamma, costs_buy=costs_buy)

        costs_sell = pd.Series(
            index=mu.index[::-1], data=np.linspace(0.0, 0.1, mu.size)
        )
        with self.assertRaises(ValueError):
            _ = mvp.efficient_portfolio(gamma, costs_sell=costs_sell)

        initial_holdings = pd.Series(
            index=mu.index[::-1], data=np.linspace(0.0, 0.01, mu.size)
        )
        with self.assertRaises(ValueError):
            _ = mvp.efficient_portfolio(gamma, initial_holdings=initial_holdings)

    def test_risk_free_asset(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 12.5

        # Parameters chosen such that some fraction is invested in cash
        mvp = MeanVariancePortfolio(mu, cov_matrix)
        pf = mvp.efficient_portfolio(gamma, rf_return=0.0025)

        self.assertTrue(hasattr(pf, "x_rf"))
        self.assertGreater(pf.x_rf, 0.1)
        self.assertAlmostEqual(pf.ret, mu @ pf.x + 0.0025 * pf.x_rf)
