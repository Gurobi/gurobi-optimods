import gurobipy as gp
from gurobipy import GRB
import numpy as np
from dataclasses import dataclass


@dataclass
class QuboResult:
    solution: np.ndarray
    objective_value: float


def callback(model, where):

    if where == GRB.Callback.MIP:
        runtime = model.cbGet(GRB.Callback.RUNTIME)
        if runtime >= model._next_output_time:
            primal_bound = model.cbGet(GRB.Callback.MIP_OBJBST)
            dual_bound = model.cbGet(GRB.Callback.MIP_OBJBND)
            print(
                f"Time: {runtime:.0f}s, "
                f"best objective: {primal_bound:.2f}, "
                f"best bound: {dual_bound:.2f}, "
                f"gap: {100.0*(primal_bound - dual_bound)/abs(primal_bound):.2f}% "
                f"(use Ctrl+C to interrupt)"
            )
            model._next_output_time += 5

    elif where == GRB.Callback.MIPSOL:
        obj = model.cbGet(GRB.Callback.MIPSOL_OBJ)
        print(f"New QUBO solution found with objective {obj}")


def solve_qubo(
    coeff_matrix, time_limit=GRB.INFINITY, output=False, log_file=""
) -> QuboResult:
    """
    Solve a quadratic unconstrained binary optimization (QUBO) problem,
    i.e., minimize quadratic function :math:`x'Qx` defined by coefficient matrix :math:`Q`
    over a binary decision variable vector :math:`x`

    :param coeffMatrix: Quadratic coefficient matrix
    :type coeffMatrix: :class:`numpy.ndarray` or :class:`scipy.sparse`
    :param timeLimit: Time limit in seconds
    :type timeLimit: :class:`int`
    :param output: Enable progress output
    :type output: :class:`bool`
    :param logFile: Filename for Gurobi log output
    :type logFile: :class:`string`
    :return: 0/1 solution array, objective value
    :rtype: :class:`QuboResult`
    """

    if coeff_matrix.ndim != 2:
        raise ValueError("Matrix is not 2-dimensional.")

    shape = coeff_matrix.shape
    if shape[0] != shape[1]:
        raise ValueError("Matrix is not quadratic.")

    n = shape[0]

    params = {"TimeLimit": time_limit, "LogToConsole": 0, "LogFile": log_file}

    with gp.Env(params=params) as env, gp.Model(env=env) as model:

        x = model.addMVar(n, vtype=GRB.BINARY)
        model.setObjective(x @ coeff_matrix @ x, GRB.MINIMIZE)

        if output:
            model._next_output_time = 5
            model.optimize(callback)
        else:
            model.optimize()

        if model.SolCount == 0:
            raise ValueError(
                "No solution found, potentially because of a very low time limit."
            )

        return QuboResult(solution=x.X.round(), objective_value=model.ObjVal)
