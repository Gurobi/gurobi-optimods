import unittest


class TestExamples(unittest.TestCase):
    @unittest.expectedFailure  # not implemented yet
    def test_card_regression(self):
        import examples.card_regression

    def test_workforce(self):
        import examples.workforce

    def test_bipartite_matching(self):
        import examples.bipartite_matching

    def test_weighted_matching(self):

        import examples.weighted_matching
