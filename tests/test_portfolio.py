import unittest
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

from gurobi_optimods.datasets import load_portfolio
from gurobi_optimods.portfolio import MeanVariancePortfolio


class TestMeanVariancePortfolio(unittest.TestCase):
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
        x = mvp.efficient_portfolio(0.5)

    def test_two_assets(self):
        cov_matrix = np.array([[3, 0.5], [0.5, 2]])
        mu = np.array([1, -0.1])
        mvp = MeanVariancePortfolio(mu, cov_matrix)
        x = mvp.efficient_portfolio(0.5)
        assert_allclose(x, [0.925, 0.075], atol=1e-6)

    def test_number_of_trades(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0

        mvp = MeanVariancePortfolio(mu, cov_matrix)
        x_unconstrained = mvp.efficient_portfolio(gamma)
        self.assertGreater((x_unconstrained > 1e-4).sum(), 3)

        x_3 = mvp.efficient_portfolio(gamma, max_trades=3)
        self.assertLessEqual((x_3 > 1e-4).sum(), 3)

    def test_transaction_fees_long(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0
        fees_buy = 1e-4

        mvp = MeanVariancePortfolio(mu, cov_matrix)
        x = mvp.efficient_portfolio(gamma, fees_buy=fees_buy)
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
        x = mvp.efficient_portfolio(gamma, max_total_short=0.3)
        n_trades = (x < 0).sum()
        self.assertGreater(n_trades, 0)
        self.assertAlmostEqual(x.sum(), 1)

        # Ensure that transaction fees are paid out of the portfolio
        x = mvp.efficient_portfolio(gamma, max_total_short=0.3, fees_sell=fees_sell)
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
        x = mvp.efficient_portfolio(gamma, min_long=0.0)
        small_trades = (x > 1e-6) & (x < (0.03 - 1e-6))
        self.assertGreater(small_trades.sum(), 0)

        # Ensure that we _don't_ have a small trade w minimum buy in
        x = mvp.efficient_portfolio(gamma, min_long=0.03)
        small_trades = (x > 1e-6) & (x < (0.03 - 1e-6))
        self.assertEqual(small_trades.sum(), 0)

    def test_max_total_short(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0

        mvp = MeanVariancePortfolio(mu, cov_matrix)

        # Ensure that by default we don't go short
        x = mvp.efficient_portfolio(gamma)
        self.assertEqual((x < 0).sum(), 0)

        # Ensure that we take advantage of leverage
        x = mvp.efficient_portfolio(gamma, max_total_short=0.1)
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
        x = mvp.efficient_portfolio(gamma)
        self.assertEqual((x < 0).sum(), 0)

        # Adding a min_short constraint doesn't change a thing
        x_other = mvp.efficient_portfolio(gamma, min_short=0.01)
        assert_allclose(x.to_numpy(), x_other.to_numpy())

    def test_min_short_1(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0

        mvp = MeanVariancePortfolio(mu, cov_matrix)
        x = mvp.efficient_portfolio(gamma, max_total_short=0.5)

        # Ensure that we do have tiny short positions
        small_trades = (x < -1e-6) & (x > (-0.05 + 1e-6))
        self.assertGreater(small_trades.sum(), 0)

        # Ensure that we still have short positions but beyond the threshold
        x = mvp.efficient_portfolio(gamma, max_total_short=0.5, min_short=0.05)
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
        x_with = mvp.efficient_portfolio(gamma, initial_holdings=x0)
        x_without = mvp.efficient_portfolio(gamma, initial_holdings=None)
        assert_array_equal(x_with, x_without)

    def test_start_portfolio_without_restrictions(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0

        # If there are no additional restrictions, the resulting portfolio
        # should be the same as without initial holdings sunk-cost-fallacy
        mvp = MeanVariancePortfolio(mu, cov_matrix)
        x0 = 1.0 / mu.size * np.ones(mu.size)
        x_with = mvp.efficient_portfolio(gamma, initial_holdings=x0)
        x_without = mvp.efficient_portfolio(gamma, initial_holdings=None)
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
        x_with = mvp.efficient_portfolio(gamma, initial_holdings=x0)
        x_without = mvp.efficient_portfolio(gamma, initial_holdings=None)
        assert_allclose(x_with, x_without, atol=1e-6)

    def test_start_portfolio_limit_trades(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0

        mvp = MeanVariancePortfolio(mu, cov_matrix)
        x0 = 1.0 / mu.size * np.ones(mu.size)
        # Ensure that we use more than 3 trades if there is no limit
        x1_nomaxtrades = mvp.efficient_portfolio(gamma, initial_holdings=x0)
        trades = ((x1_nomaxtrades - x0) > 1e-4) | ((x1_nomaxtrades - x0) < -1e-4)
        self.assertGreater(trades.sum(), 3)

        # Ensure that we only use 3 trades with the limit
        x1_maxtrades = mvp.efficient_portfolio(gamma, initial_holdings=x0, max_trades=3)
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
        )

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
        # Ensure that we _do_ have a small trade w/o minimum buy
        x = mvp.efficient_portfolio(gamma, min_long=0.0, initial_holdings=x0)
        trades = x - x0
        small_trades = (trades > 1e-6) & (trades < (0.03 - 1e-6))
        self.assertGreater(small_trades.sum(), 0)

        # Ensure that we _don't_ have a small trade w/ minimum buy
        x = mvp.efficient_portfolio(gamma, min_long=0.03, initial_holdings=x0)
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
        x0 = mvp.efficient_portfolio(gamma, max_total_short=0.1)

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
        )
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
        )
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
        )
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
        x = mvp.efficient_portfolio(gamma, max_positions=3)
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
        )

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

        x = mvp.efficient_portfolio(
            gamma,
            initial_holdings=x0,
            max_positions=6,
            max_trades=3,
        )
        # This needs to be infeasible.
        self.assertIsNone(x)

    def test_transaction_costs_long(self):
        data = load_portfolio()
        cov_matrix = data.cov()
        mu = data.mean()
        gamma = 100.0
        costs = 0.0025

        mvp = MeanVariancePortfolio(mu, cov_matrix)
        x = mvp.efficient_portfolio(gamma, costs_buy=costs)
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
        x = mvp.efficient_portfolio(gamma, costs_buy=costs, fees_buy=fees)
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
        )
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
        )
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
        )
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
        )

        self.assertAlmostEqual(x[x < 0].sum(), -0.1)
