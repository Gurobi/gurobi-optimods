# Implementation of your new mod. This should be copied to
# src/gurobi_optimods/<mod-name>.py. You may alternatively want to include
# your mod in an existing file, if it coexists naturally with other mods.
#
# In general the public API should be a single class or function. Go with
# whatever makes the most sense for this mod.

import gurobipy as gp
from gurobipy import GRB
import numpy as np


def solve_qubo(coeff_matrix):
    """
    A sphinx-compatible docstring

    :param data1: Data structure for first argument
    :type data1: pd.DataFrame
    """

    if coeff_matrix is None:
        return None

    # TODO return useful error
    if coeff_matrix.ndim != 2:
        return None

    shape = coeff_matrix.shape
    if shape[0] != shape[1]:
        return None

    n = shape[0]

    with gp.Env() as env, gp.Model(env=env) as model:

        x = model.addMVar(n, vtype=GRB.BINARY)
        model.setObjective(x @ coeff_matrix @ x, GRB.MINIMIZE)

        model.optimize()

        return x.X.round()
