"""
(Elementary) Shortest Path Problem with Resource Constraints
------------------------------------------------------------
"""

import logging
from dataclasses import dataclass

import gurobipy as gp
import numpy as np
from gurobipy import GRB

from gurobi_optimods.utils import optimod

logger = logging.getLogger(__name__)


@dataclass
class ShortestPath:
    """
    Solution to a QUBO problem.

    Attributes
    ----------
    solution : ndarray
        0/1 array of variable values in the solution
    objective_value : float
        The objective function value for this solution
    """

    nodes: list
    cost: float


@optimod()
def shortest_path(
    node_data, arc_data, source, target, limits, elementary=True, *, create_env
) -> ShortestPath:
    """
    An optimod that solves an important problem

    :param data: Description of argument
    :type data: Type of argument

    ... describe additional arguments ...

    :return: Description of returned result
    :rtype: Type of returned result
    """

    # if coeff_matrix.ndim != 2:
    #     raise ValueError("Matrix is not 2-dimensional.")

    params = {"LogToConsole": 0}

    with create_env(params=params) as env, gp.Model(env=env) as model:
        model.optimize()

        if model.SolCount == 0:
            raise ValueError(
                "No solution found, potentially because of a very low time limit."
            )

        return ShortestPath(nodes=[0], cost=model.ObjVal)
