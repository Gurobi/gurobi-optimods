"""
(Elementary) Shortest Path Problem with Resource Constraints
------------------------------------------------------------
"""

import logging
from dataclasses import dataclass

import gurobipy as gp
import pandas as pd
from gurobipy import GRB

from gurobi_optimods.utils import optimod

logger = logging.getLogger(__name__)


@dataclass
class ShortestPath:
    """
    Solution to a QUBO problem.

    Attributes
    ----------
    nodes : list
        ordered list of visited nodes along the path
    cost : float
        total path costs
    """

    nodes: list
    cost: float


@optimod()
def shortest_path(
    node_data, arc_data, source, target, limits={}, elementary=True, *, create_env
) -> ShortestPath:
    """
    An optimod that solves an important problem

    :param data: Description of argument
    :type data: Type of argument

    ... describe additional arguments ...

    :return: Description of returned result
    :rtype: Type of returned result
    """

    # add node costs to arc costs

    params = {"LogToConsole": 0}

    with create_env(params=params) as env, gp.Model(env=env) as model:
        model.ModelSense = GRB.MINIMIZE

        arc_df = arc_data.gppd.add_vars(
            model, vtype=GRB.INTEGER, lb=0, obj="cost", name="arc"
        )

        node_flows = pd.DataFrame(
            {
                "inflow": arc_df["flow"].groupby("to").sum(),
                "outflow": arc_df["flow"].groupby("from").sum(),
            }
        ).fillna(0)

        # start path at source, end path at target

        # flow conservation constraints (TODO exclude source and target)
        node_flows.gppd.add_constrs(model, "inflow == outflow", name="balance")

        # restriction to elementary paths

        # resource restrictions

        logger.info(f"Solving shortest path problem")

        model.optimize()

        if model.Status in [GRB.INFEASIBLE, GRB.INF_OR_UNBD]:
            raise ValueError("Infeasible or unbounded")

        if model.SolCount == 0:
            raise ValueError(
                "No solution found, potentially because of a very low time limit."
            )

        return ShortestPath(nodes=[0], cost=model.ObjVal)
