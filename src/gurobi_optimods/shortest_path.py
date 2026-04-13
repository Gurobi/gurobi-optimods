"""
(Elementary) Shortest Path Problem with Resource Constraints
------------------------------------------------------------
"""

import logging
from dataclasses import dataclass

import gurobipy as gp
import gurobipy_pandas as gppd
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

    # add node properties to arc properties (matching the source node of each arc)
    print(node_data)
    print(arc_data)

    arcnode_data = arc_data.add(node_data, fill_value=0)
    print(arcnode_data)

    params = {"LogToConsole": 1}

    with create_env(params=params) as env, gp.Model(env=env) as model:
        model.ModelSense = GRB.MINIMIZE

        arc_df = arcnode_data.gppd.add_vars(
            model, vtype=GRB.INTEGER, lb=0, obj="cost", name="arc"
        )

        # add cost of target node
        model.ObjCon = node_data.loc[target, "cost"]

        node_flows = pd.DataFrame(
            {
                "inflow": arc_df["arc"].groupby("j").sum(),
                "outflow": arc_df["arc"].groupby("i").sum(),
            }
        ).fillna(0)
        # balance value: source = 1, target = -1, 0 otherwise
        node_flows["balance"] = (node_flows.index == source).astype(int) - (
            node_flows.index == target
        ).astype(int)

        model.update()
        print(node_flows)

        # flow conservation constraints
        node_flows.gppd.add_constrs(
            model, "outflow - inflow == balance", name="balance_constr"
        )

        # restriction to elementary paths
        if elementary:
            node_flows.gppd.add_constrs(model, "inflow <= 1", name="elementary_constr")

        # TODO resource restrictions

        logger.info(f"Solving shortest path problem")

        model.write("model.lp")
        model.optimize()

        if model.Status in [GRB.INFEASIBLE, GRB.INF_OR_UNBD]:
            raise ValueError("Infeasible or unbounded")

        if model.SolCount == 0:
            raise ValueError(
                "No solution found, potentially because of a very low time limit."
            )

        # TODO extract solution and create node visiting sequence
        arc_df["solution"] = arc_df["arc"].gppd.X
        print(arc_df)

        return ShortestPath(nodes=[0], cost=model.ObjVal)
