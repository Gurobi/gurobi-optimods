# Unit tests for your mod. This is separate from the example tests, it should
# test smaller pieces of functionality, and can also test multiple cases.
# This template should be copied to tests/test_<mod-name>.py

import unittest
import numpy as np
from numpy.testing import assert_allclose

from gurobi_optimods.datasets import load_portfolio
from gurobi_optimods.portfolio import MeanVariancePortfolio


class TestMod(unittest.TestCase):
    def test_datasets(self):
        data = load_portfolio()
        self.assertEqual(
            set(data.keys()),
            {"AA", "BB", "CC", "DD", "EE", "FF", "GG", "HH", "II", "JJ"},
        )

    def test_example_data(self):
        data = load_portfolio()
        Sigma = data.cov()
        mu = data.mean()

        mvp = MeanVariancePortfolio(Sigma, mu)
        x = mvp.efficient_portfolio(0.5)

    def test_two_assets(self):
        Sigma = np.array([[3, 0.5], [0.5, 2]])
        mu = np.array([1, -0.1])
        mvp = MeanVariancePortfolio(Sigma, mu)
        x = mvp.efficient_portfolio(0.5)
        assert_allclose(x, [0.925, 0.075], atol=1e-6)

    def test_number_of_trades(self):
        data = load_portfolio()
        Sigma = data.cov()
        mu = data.mean()
        gamma = 100.0

        mvp = MeanVariancePortfolio(Sigma, mu)
        x_unconstrained = mvp.efficient_portfolio(gamma)
        self.assertGreater((x_unconstrained > 1e-4).sum(), 3)

        x_3 = mvp.efficient_portfolio(gamma, max_trades=3)
        self.assertLessEqual((x_3 > 1e-4).sum(), 3)

    def test_transaction_fees(self):
        data = load_portfolio()
        Sigma = data.cov()
        mu = data.mean()
        gamma = 100.0
        fees_buy = 1e-4

        mvp = MeanVariancePortfolio(Sigma, mu)
        x = mvp.efficient_portfolio(gamma, fees_buy=fees_buy)
        n_trades = (x > 1e-4).sum()
        self.assertLessEqual(x.sum(), 1 - n_trades * fees_buy)

    def test_minimum_buy_in(self):
        data = load_portfolio()
        Sigma = data.cov()
        mu = data.mean()
        gamma = 100.0

        mvp = MeanVariancePortfolio(Sigma, mu)

        # Ensure that we _have_ a small trade w/o minimum buy in
        x = mvp.efficient_portfolio(gamma, min_buy_in=0.0)
        small_trades = (x > 1e-6) & (x < (0.03 - 1e-6))
        self.assertGreater(small_trades.sum(), 0)

        # Ensure that we _don't_ have a small trade w minimum buy in
        x = mvp.efficient_portfolio(gamma, min_buy_in=0.03)
        small_trades = (x > 1e-6) & (x < (0.03 - 1e-6))
        self.assertEqual(small_trades.sum(), 0)
