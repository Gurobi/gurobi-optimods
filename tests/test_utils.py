import io
import os
import tempfile
import unittest
import warnings
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

import gurobipy as gp

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

        buffer_stdout = io.StringIO()
        buffer_stderr = io.StringIO()

        with redirect_stdout(buffer_stdout), redirect_stderr(buffer_stderr):
            self.mod()

        buffer_stdout.seek(0)
        buffer_stderr.seek(0)

        self.assertIn("Gurobi Optimizer", buffer_stdout.read())
        self.assertEqual(buffer_stderr.read(), "")

    def test_silent(self):
        # silent=True disables all output

        buffer_stdout = io.StringIO()
        buffer_stderr = io.StringIO()

        with redirect_stdout(buffer_stdout), redirect_stderr(buffer_stderr):
            self.mod(silent=True)

        buffer_stdout.seek(0)
        buffer_stderr.seek(0)

        self.assertEqual(buffer_stdout.read(), "")
        self.assertEqual(buffer_stderr.read(), "")

    def test_logfile(self):
        # Write to a target log file

        buffer_stdout = io.StringIO()
        buffer_stderr = io.StringIO()

        with tempfile.TemporaryDirectory() as tempdir, redirect_stdout(
            buffer_stdout
        ), redirect_stderr(buffer_stderr):

            logfile = os.path.join(tempdir, "tmp.log")
            self.mod(logfile=logfile)

            logfile_text = Path(logfile).read_text()

        buffer_stdout.seek(0)
        buffer_stderr.seek(0)

        self.assertIn("Gurobi Optimizer", buffer_stdout.read())
        self.assertEqual(buffer_stderr.read(), "")
        self.assertIn("Gurobi Optimizer", logfile_text)

    def test_logfile_closed(self):
        # Ensure no resource warnings due to files left open

        with warnings.catch_warnings(record=True, category=ResourceWarning) as w:
            self.mod(logfile="mod.log")
            assert not w
