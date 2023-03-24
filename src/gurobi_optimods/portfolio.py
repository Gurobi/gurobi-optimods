import gurobipy as gp

# In general the public API should be a single class or function. Go with
# whatever makes the most sense for this mod.

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
