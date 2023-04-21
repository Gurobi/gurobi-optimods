import gurobipy as gp
from gurobipy import GRB
import numpy as np
from dataclasses import dataclass


@dataclass
class QuboResult:
    solution: np.array
    objective_value: float


def solve_qubo(coeff_matrix) -> QuboResult:
    """
    Solve a quadratic unconstrained binary optimization (QUBO) problem,
    i.e., minimize quadratic function x'Qx defined by coefficient matrix Q
    over a binary decision variable vector x

    :param coeff_matrix: quadratic coefficient matrix
    :type coeff_matrix: np.array or scipy.sparse
    """

    if coeff_matrix is None:
        return None

    if coeff_matrix.ndim != 2:
        raise ValueError("Matrix is not 2-dimensional.")

    shape = coeff_matrix.shape
    if shape[0] != shape[1]:
        raise ValueError("Matrix is not quadratic.")

    n = shape[0]

    with gp.Env() as env, gp.Model(env=env) as model:

        x = model.addMVar(n, vtype=GRB.BINARY)
        model.setObjective(x @ coeff_matrix @ x, GRB.MINIMIZE)

        model.optimize()

        return QuboResult(objective_value=model.ObjVal, solution=x.X.round())
