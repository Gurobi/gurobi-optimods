import unittest

import numpy as np
from numpy.testing import assert_allclose

from gurobi_optimods.regression import LADRegression


class TestLADRegression(unittest.TestCase):
    def test_two_points(self):
        # L1 regression through two points in 2D creates a perfect fit
        X_train = np.array([[1.0], [3.0]])
        y_train = np.array([2.0, 6.0])
        reg = LADRegression()
        reg.fit(X_train, y_train)

        # Known fit
        assert_allclose(reg.intercept_, 0.0)
        assert_allclose(reg.coef_, np.array([2.0]))

        # Test predictions
        y_pred = reg.predict(np.array([[1.0], [2.0], [3.0]]))
        assert_allclose(y_pred, np.array([2.0, 4.0, 6.0]))

    def test_random(self):
        # Plug in some random data and check shape consistency
        X_train = np.random.random((100, 5))
        y_train = np.random.random((100))
        reg = LADRegression()
        reg.fit(X_train, y_train)

        # Check predictions are the right shape
        y_pred = reg.predict(np.random.random((30, 5)))
        self.assertEqual(y_pred.shape, (30,))
