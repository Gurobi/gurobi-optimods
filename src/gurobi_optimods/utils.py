# One possible idea for handling output suppression/file logging. Mods need
# to be decorated with @optimod and accept a create_env keyword arguemnt. This
# factory function should be used to create the mod's gurobi environments.
#
#   @optimod(mod_logger=logger)  # pass the mod's module logger
#   def my_mod(inputs, *, create_env):
#       with create_env() as env, gp.Model(env=env) as model:
#           # do stuff with model
#
# The mod then gets two arguments for free: `silent`, where `silent=True`
# suppresses all output, and logfile=<file-path> creates a log file including
# logger output from the mod and gurobi output.
#
# FIXME this won't handle multi-threaded stuff well (context manipulates
# global logger handlers). But it's simpler to implement in the mods than
# passing some funky log collecting object around.

import functools
import logging
import re
from contextlib import contextmanager
from typing import Optional

import gurobipy as gp

global_mod_logger = logging.getLogger(r"gurobi_optimods")
grb_logger = logging.getLogger(r"gurobipy")
re_module_base_name = re.compile(r"gurobipy\.|gurobi_optimods\.")


class ShortFormatter(logging.Formatter):
    def __init__(self):
        super().__init__("%(name)s: %(message)s")

    def format(self, record: logging.LogRecord):
        record.name = re_module_base_name.sub("", record.name, 1)
        s = super().format(record)
        return s


@contextmanager
def _mod_context(
    *, mod_logger: logging.Logger, log_to_console: bool, log_to_file: Optional[str]
):
    if not log_to_console and log_to_file:
        raise ValueError("Cannot disable console output and log to file")

    # Base setting: silence
    grbenv_params = {"OutputFlag": 0}

    if log_to_console:
        # Gurobi console output handled by environment
        grbenv_params["OutputFlag"] = 1

        # Send mod logs to console
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(logging.Formatter("%(message)s"))
        mod_logger.setLevel(logging.INFO)
        mod_logger.addHandler(ch)

    if log_to_file:
        # Handle all file logging using the python logger
        assert grbenv_params["OutputFlag"] == 1

        fh = logging.FileHandler(log_to_file)
        fh.setLevel(logging.INFO)
        fh.setFormatter(ShortFormatter())

        # Send both mod logs and gurobi logs to the same file handler
        mod_logger.addHandler(fh)
        grb_logger.setLevel(logging.INFO)
        grb_logger.addHandler(fh)

    # Environment factory for decorated mod to use
    def create_env():
        return gp.Env(params=grbenv_params)

    try:
        yield create_env

    finally:
        if log_to_console:
            mod_logger.removeHandler(ch)

        if log_to_file:
            mod_logger.removeHandler(fh)
            grb_logger.removeHandler(fh)
            fh.close()


def optimod(mod_logger=None):
    if mod_logger is None:
        mod_logger = global_mod_logger

    def optimod_decorator(func):
        @functools.wraps(func)
        def optimod_decorated(*args, silent=False, logfile=None, **kwargs):
            with _mod_context(
                mod_logger=mod_logger,
                log_to_console=not silent,
                log_to_file=logfile,
            ) as create_env:
                return func(*args, create_env=create_env, **kwargs)

        return optimod_decorated

    return optimod_decorator
