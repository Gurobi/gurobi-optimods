import unittest

from numpy.testing import assert_allclose, assert_array_equal
from pandas.testing import assert_frame_equal


class TestL1Regression(unittest.TestCase):
    def test_compare_examples(self):

        import examples.l1_regression.gurobipy as gp_impl
        import examples.l1_regression.nupstup as ns_impl

        assert_allclose(gp_impl.y_pred, ns_impl.y_pred)


class TestWorkforce(unittest.TestCase):
    def test_compare_examples(self):

        import examples.workforce.gurobipy as gp_impl
        import examples.workforce.nupstup as ns_impl

        assert_frame_equal(
            gp_impl.assigned_shifts.reset_index(drop=True), ns_impl.assigned_shifts
        )


class TestBipartiteMatching(unittest.TestCase):
    def test_compare_examples(self):

        import examples.bipartite_matching.gurobipy as gp_impl
        import examples.bipartite_matching.nupstup as ns_impl

        assert_array_equal(gp_impl.matching.todense(), ns_impl.matching.todense())


class TestWeightedMatching(unittest.TestCase):
    def test_compare_examples(self):

        import examples.weighted_matching.gurobipy as gp_impl
        import examples.weighted_matching.nupstup as ns_impl

        print(gp_impl.matching.toarray())
        print(ns_impl.matching.toarray())

        assert_array_equal(gp_impl.matching.todense(), ns_impl.matching.todense())
