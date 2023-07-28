"""
Mods must be decorated with @optimod and accept a create_env keyword argument.
This factory function should be used to create the mod's gurobi environments::

   @optimod()
   my_mod(inputs, *, create_env):
       with create_env() as env, gp.Model(env=env) as model:
           # formulate and solve model

The mod then gets some arguments for free:

- `verbose`, where `verbose=False` suppresses all output;
- logfile=<file-path>, which writes output to a log file; and
- solver_params, which the Mod caller can use to pass parameters to Gurobi.

Parameters can also be passed as a dictionary to create_env if the Mod requires
some specific settings.

Note that this captures output via the gurobipy and optimod python loggers. It
may not work as expected when multithreading in Python.
"""

import functools
import logging
import re
import sys
from contextlib import contextmanager
from typing import Dict, Optional

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
    *,
    mod_logger: logging.Logger,
    log_to_console: bool,
    log_to_file: Optional[str],
    user_params: Optional[Dict],
):
    if not log_to_console and log_to_file:
        raise ValueError("Cannot disable console output and log to file")

    # Base setting: silence
    decorator_params = {"OutputFlag": 0}

    if log_to_console:
        # Gurobi console output handled by environment
        decorator_params["OutputFlag"] = 1

        # Send mod logs to console
        ch = logging.StreamHandler(stream=sys.stdout)
        ch.setLevel(logging.INFO)
        ch.setFormatter(logging.Formatter("%(message)s"))
        mod_logger.setLevel(logging.INFO)
        mod_logger.addHandler(ch)

    if log_to_file:
        # Handle all file logging using the python logger
        assert decorator_params["OutputFlag"] == 1

        fh = logging.FileHandler(log_to_file)
        fh.setLevel(logging.INFO)
        fh.setFormatter(ShortFormatter())

        # Send both mod logs and gurobi logs to the same file handler
        mod_logger.addHandler(fh)
        grb_logger.setLevel(logging.INFO)
        grb_logger.addHandler(fh)

    # Environment factory for decorated mod to use
    def create_env(params=None):
        final_params = {}
        final_params.update(decorator_params)
        if params:
            final_params.update(params)
        if user_params:
            final_params.update(user_params)
        return gp.Env(params=final_params)

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
        def optimod_decorated(
            *args, verbose=True, logfile=None, solver_params=None, **kwargs
        ):
            with _mod_context(
                mod_logger=mod_logger,
                log_to_console=verbose,
                log_to_file=logfile,
                user_params=solver_params,
            ) as create_env:
                try:
                    return func(*args, create_env=create_env, **kwargs)

                except gp.GurobiError as ge:
                    if ge.errno == gp.GRB.ERROR_SIZE_LIMIT_EXCEEDED:
                        pass  # fall through
                    else:
                        raise

                # We can only fall through to here due to SIZE_LIMIT_EXCEEDED,
                # so raise a more optimods-appropriate error. Raise here
                # (instead of directly in the except block above) to avoid a
                # confusing double stack trace.
                raise ValueError(
                    "Given data exceeds Gurobi trial license limits; please see "
                    "https://support.gurobi.com/hc/en-us/articles/15801588452241 "
                    "to resolve this issue"
                )

        optimod_decorated._decorated_mod = True
        return optimod_decorated

    return optimod_decorator


def run(f):
    """A decorator that will run the body of the decorated function,
    and then store the result of that function in place of the function’s name.

    From: https://chrisjrn.com/2021/09/10/talk-notes--on-the-use-and-misuse-of-decorators/#run-the-scoped-block
    """
    return f()
