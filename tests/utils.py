"""
Testing utilities, not to be confused with unit tests of gurobi_optimods.utils
found in tests/test_utils.py
"""

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
