import unittest
import numpy as np
import pandas as pd

from gurobi_optimods.datasets import load_sharpe_ratio
from gurobi_optimods.sharpe_ratio import max_sharpe_ratio


class TestSharpeRatio(unittest.TestCase):
    def test_dataset_keys(self):
        data = load_sharpe_ratio()
        assets = {
            "A",
            "B",
            "C",
            "D",
            "E",
            "F",
            "G",
            "H",
            "I",
            "J",
            "K",
            "L",
            "M",
            "N",
            "O",
            "P",
            "Q",
            "R",
            "S",
            "T",
            "U",
            "V",
            "W",
            "X",
            "Y",
            "Z",
        }
        self.assertEqual(set(data.cov_matrix.keys()), assets)
        self.assertEqual(set(data.mu.keys()), assets)

    def test_invalid_arg_types(self):
        data = load_sharpe_ratio()

        with self.assertRaises(TypeError):
            max_sharpe_ratio(None, None, None)

        with self.assertRaises(TypeError):
            max_sharpe_ratio(1, data.mu)

        with self.assertRaises(TypeError):
            max_sharpe_ratio(data.mu, data.mu)

        with self.assertRaises(TypeError):
            max_sharpe_ratio(data.cov_matrix, data.cov_matrix)

        with self.assertRaises(TypeError):
            max_sharpe_ratio(data.cov_matrix, 1)

        with self.assertRaises(TypeError):
            max_sharpe_ratio(data.cov_matrix, data.mu, "0")

    def test_incorrect_numpy_dimensions(self):
        with self.assertRaises(ValueError):
            max_sharpe_ratio(np.ones((1, 1, 1)), np.ones(1))

        with self.assertRaises(ValueError):
            max_sharpe_ratio(np.ones(1), np.ones(1))

        with self.assertRaises(ValueError):
            max_sharpe_ratio(np.ones((1, 1)), np.ones((1, 1)))

    def test_numpy_inputs(self):
        data = load_sharpe_ratio()

        portfolio, ratio = max_sharpe_ratio(
            data.cov_matrix.to_numpy(), data.mu.to_numpy()
        )
        self.assertIsInstance(portfolio, np.ndarray)
        self.assertEqual(portfolio.shape, data.mu.to_numpy().shape)
        self.assertAlmostEqual(portfolio.sum(), 1, delta=1e-6)
        self.assertIsInstance(ratio, float)

    def test_pandas_inputs(self):
        data = load_sharpe_ratio()

        portfolio, ratio = max_sharpe_ratio(data.cov_matrix, data.mu)
        self.assertIsInstance(portfolio, pd.Series)
        self.assertTrue(data.cov_matrix.index.identical(portfolio.index))
        self.assertAlmostEqual(portfolio.sum(), 1, delta=1e-6)
        self.assertIsInstance(ratio, float)

        portfolio, ratio = max_sharpe_ratio(data.cov_matrix, data.mu.to_numpy())
        self.assertIsInstance(portfolio, pd.Series)
        self.assertTrue(data.cov_matrix.index.identical(portfolio.index))
        self.assertAlmostEqual(portfolio.sum(), 1, delta=1e-6)
        self.assertIsInstance(ratio, float)

        portfolio, ratio = max_sharpe_ratio(data.cov_matrix.to_numpy(), data.mu)
        self.assertIsInstance(portfolio, pd.Series)
        self.assertTrue(data.cov_matrix.index.identical(portfolio.index))
        self.assertAlmostEqual(portfolio.sum(), 1, delta=1e-6)
        self.assertIsInstance(ratio, float)

    def test_mismatched_indices(self):
        data = load_sharpe_ratio()

        mu = data.mu.rename({"A": "foo"})
        with self.assertRaises(ValueError):
            max_sharpe_ratio(data.cov_matrix, mu)

        mu = data.mu.reindex(data.mu.index[::-1])
        with self.assertRaises(ValueError):
            max_sharpe_ratio(data.cov_matrix, mu)

    def test_single_asset(self):
        cov_matrix = np.eye(1)
        mu = np.array([0.1])
        portfolio, ratio = max_sharpe_ratio(cov_matrix, mu)
        self.assertIsInstance(portfolio, np.ndarray)
        self.assertEqual(portfolio.shape, (1,))
        self.assertAlmostEqual(portfolio.sum(), 1, delta=1e-6)
        self.assertIsInstance(ratio, float)
        self.assertAlmostEqual(ratio, 0.1 / cov_matrix[0][0] ** 0.5, delta=1e-6)

    def test_rf_rate(self):
        cov_matrix = np.eye(1)
        mu = np.array([0.1])
        rf_rate = 0.01
        portfolio, ratio = max_sharpe_ratio(cov_matrix, mu, rf_rate)
        self.assertIsInstance(portfolio, np.ndarray)
        self.assertEqual(portfolio.shape, (1,))
        self.assertAlmostEqual(portfolio.sum(), 1, delta=1e-6)
        self.assertIsInstance(ratio, float)
        self.assertAlmostEqual(
            ratio, (mu[0] - rf_rate) / cov_matrix[0][0] ** 0.5, delta=1e-6
        )

    def test_risk_free_investment_is_best(self):
        cov_matrix = np.eye(3)
        mu = np.array([-0.1, -0.01, -0.3])

        with self.assertRaises(ValueError):
            max_sharpe_ratio(cov_matrix, mu)

        mu = np.array([0.04, 0.03, 0.02])
        rf_rate = 0.041
        with self.assertRaises(ValueError):
            max_sharpe_ratio(cov_matrix, mu, rf_rate)

    def test_single_asset_with_positive_return(self):
        cov_matrix = np.eye(3)
        mu = np.array([-0.1, -0.01, 0.01])
        portfolio, ratio = max_sharpe_ratio(cov_matrix, mu)
        self.assertIsInstance(portfolio, np.ndarray)
        self.assertEqual(portfolio.shape, (3,))
        self.assertAlmostEqual(portfolio.sum(), 1, delta=1e-6)
        self.assertIsInstance(ratio, float)
        self.assertAlmostEqual(ratio, mu[2] / cov_matrix[0][0] ** 0.5, delta=1e-6)

    def test_dataset_maximal_ratio(self):
        data = load_sharpe_ratio()

        portfolio, ratio = max_sharpe_ratio(data.cov_matrix, data.mu)
        self.assertIsInstance(portfolio, pd.Series)
        self.assertTrue(data.cov_matrix.index.identical(portfolio.index))
        self.assertAlmostEqual(portfolio.sum(), 1, delta=1e-6)
        self.assertIsInstance(ratio, float)
        self.assertAlmostEqual(ratio, 2.417932187428722, delta=1e-6)
