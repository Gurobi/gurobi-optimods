import unittest

import numpy as np

from gurobi_optimods.regression import LADRegression


class TestLADRegression(unittest.TestCase):
    def test_random(self):
        reg = LADRegression()
        X_train = np.random.random((100, 5))
        y_train = np.random.random((100))
        reg.fit(X_train, y_train)
        y_pred = reg.predict(np.random.random((30, 5)))
        self.assertEqual(y_pred.shape, (30,))
