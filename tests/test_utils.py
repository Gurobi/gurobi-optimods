"""
Test of gurobi_optimods.utils, not to be confused with utilities for testing
found in tests/utils.py
"""

import io
import os
import tempfile
import unittest
import warnings
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

import gurobipy as gp
from gurobipy import GRB

from gurobi_optimods.utils import optimod


class TestOptimodDecorator(unittest.TestCase):
    def setUp(self):
        @optimod()
        def mod(*, create_env):
            with create_env() as env, gp.Model(env=env) as model:
                model.optimize()

        self.mod = mod

    def test_basic(self):
        # By default, logs are written to standard output

        with redirect_stdout(io.StringIO()) as buffer_stdout, redirect_stderr(
            io.StringIO()
        ) as buffer_stderr:
            self.mod()

        self.assertIn("Gurobi Optimizer", buffer_stdout.getvalue())
        self.assertEqual(buffer_stderr.getvalue(), "")

    def test_not_verbose(self):
        # verbose=False disables all output

        with redirect_stdout(io.StringIO()) as buffer_stdout, redirect_stderr(
            io.StringIO()
        ) as buffer_stderr:
            self.mod(verbose=False)

        self.assertEqual(buffer_stdout.getvalue(), "")
        self.assertEqual(buffer_stderr.getvalue(), "")

    def test_logfile(self):
        # Write to a target log file

        with tempfile.TemporaryDirectory() as tempdir, redirect_stdout(
            io.StringIO()
        ) as buffer_stdout, redirect_stderr(io.StringIO()) as buffer_stderr:
            logfile = os.path.join(tempdir, "tmp.log")
            self.mod(logfile=logfile)
            logfile_text = Path(logfile).read_text()

        self.assertIn("Gurobi Optimizer", buffer_stdout.getvalue())
        self.assertEqual(buffer_stderr.getvalue(), "")
        self.assertIn("Gurobi Optimizer", logfile_text)

    def test_logfile_closed(self):
        # Ensure no resource warnings due to files left open

        with warnings.catch_warnings(
            record=True
        ) as w, tempfile.TemporaryDirectory() as tempdir:
            logfile = os.path.join(tempdir, "tmp.log")
            self.mod(logfile=logfile)
            assert not w

    def test_license_limit(self):
        # Ensure that we point to an appropriate license limit KB

        # Step 1: Are we using the pip license at all?
        is_limited = False
        with gp.Env() as env, gp.Model(env=env) as model:
            model.addVars(2001)
            try:
                model.optimize()
            except gp.GurobiError as ge:
                if ge.errno == gp.GRB.ERROR_SIZE_LIMIT_EXCEEDED:
                    is_limited = True

        if not is_limited:
            self.skipTest("No pip license in use")

        # Step 2: Trigger a license limit error through env creation wrapper
        @optimod()
        def too_large_mod(*, create_env):
            with create_env() as env, gp.Model(env=env) as model:
                model.addVars(2001)
                model.optimize()

        # We expect a ValueError; getting the standard license limit error is a
        # fail, everything else an error
        re_expect_message = (
            "Given data exceeds Gurobi's license limits; please see "
            "https://support.gurobi.com.*to resolve this issue"
        )
        with self.assertRaisesRegex(ValueError, re_expect_message):
            try:
                too_large_mod()
            except gp.GurobiError as ge:
                self.assertNotEqual(ge.errno, gp.GRB.ERROR_SIZE_LIMIT_EXCEEDED)
                raise


class TestOverrideParams(unittest.TestCase):
    def test_mod_override_outputflag(self):
        # The mod can pass custom parameters which override those created
        # by verbose/logfile

        @optimod()
        def mod(*, create_env):
            p = {"OutputFlag": 0}
            with create_env(params=p) as env, gp.Model(env=env) as model:
                model.optimize()

        with redirect_stdout(io.StringIO()) as buffer_stdout, redirect_stderr(
            io.StringIO()
        ) as buffer_stderr:
            # Normally, output would be produced, but the mod sets
            # outputflag=0, disabling all gurobi logging
            mod(verbose=True)

        self.assertEqual(buffer_stdout.getvalue(), "")
        self.assertEqual(buffer_stderr.getvalue(), "")

    def test_user_override_worklimit(self):
        # The user can pass through parameters which take precedence

        @optimod()
        def mod(*, create_env):
            with create_env() as env, gp.Model(env=env) as model:
                model.optimize()
                assert model.Status == GRB.WORK_LIMIT

        mod(solver_params={"WorkLimit": 0.0})
