"""
Testing utilities, not to be confused with unit tests of gurobi_optimods.utils
found in tests/test_utils.py
"""

from dataclasses import fields
import functools
import unittest

import gurobipy as gp
from gurobipy import GRB


def large_model(test_item):
    """Decorator for tests which create large models (i.e. those that the pip
    demo license does not cover). If a decorated test fails due to license limits,
    it will be skipped."""

    @functools.wraps(test_item)
    def skip_wrapper(*args, **kwargs):
        try:
            test_item(*args, **kwargs)
        except gp.GurobiError as ge:
            # Skip the test if the failure was due to licensing
            if ge.errno == GRB.Error.SIZE_LIMIT_EXCEEDED:
                raise unittest.SkipTest("Size-limited Gurobi license")
            # Otherwise, let the error go through as-is
            raise

    return skip_wrapper


def assert_tuple_compatibility(test_case, result):
    """Check that a _LimitedTuple dataclass is (mostly) backward compatible w/ tuple."""
    length = len(result)
    # check length / iterability
    test_case.assertEqual(len(result), len(tuple(result)))
    # check that original attributes can be accessed by name and index
    for index, field in enumerate(fields(result)[:length]):
        with test_case.subTest(field=field.name):
            test_case.assertIs(getattr(result, field.name), result[index])
    # check that additional attributes can only be accessed by name
    for index, field in enumerate(fields(result)[length:], start=length):
        with test_case.subTest(field=field.name):
            test_case.assertTrue(hasattr(result, field.name))
            with test_case.assertRaises(IndexError):
                result[index]


def size_limited_license():
    result = False

    try:
        import gurobipy as gp
        from gurobipy import GRB

        with gp.Env(params={"OutputFlag": 0}) as env, gp.Model(env=env) as model:
            x = model.addVars(2001)
            model.optimize()
    except gp.GurobiError as e:
        if e.errno == GRB.Error.SIZE_LIMIT_EXCEEDED:
            result = True

    return result
