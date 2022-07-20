# Add a simple unittest here which imports both the gurobipy and nupstup
# versions of the example code, and verifies that their results match.
# This ensures the gurobipy code is a faithful representation of the
# nup behaviour.

import unittest

class TestNup(unittest.TestCase):

    def test_compare(self):
        from examples.nup.gurobipy import some_result as gp_impl
        from examples.nup.nupstup import some_result as ns_impl
        self.assertEqual(gp_impl, ns_impl)
