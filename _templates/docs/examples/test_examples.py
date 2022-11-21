# Add a simple unittest here which imports the example code, and verifies
# the results.

import unittest


class TestExamples(unittest.TestCase):
    def test_compare(self):
        from examples.mod import some_result
