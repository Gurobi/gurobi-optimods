import gurobipy as gp
from gurobipy import GRB
import numpy as np


def solve_qubo(coeff_matrix) -> np.array:
    """
    Solve a quadratic unconstrained binary optimization (QUBO) problem,
    i.e., minimize a quadratic function defined by a coefficient matrix
    over a binary decision variable vector

    :param coeff_matrix: quadratic coefficient matrix
    :type coeff_matrix: np.array
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
