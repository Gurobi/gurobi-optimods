# Implementation of your new mod. This should be copied to
# src/gurobi_optimods/<mod-name>.py. You may alternatively want to include
# your mod in an existing file, if it coexists naturally with other mods.
#
# In general the public API should be a single class or function. Go with
# whatever makes the most sense for this mod.

import gurobipy as gp


def solve_mod(data):
    """
    A sphinx-compatible docstring

    :param data1: Data structure for first argument
    :type data1: pd.DataFrame
    """
    with gp.Env() as env, gp.Model(env=env) as model:
        # build model
        model.optimize()
        # post-process and return solution
        return
