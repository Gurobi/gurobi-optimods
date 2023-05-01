# Implementation of your new mod. This should be copied to
# src/gurobi_optimods/<mod-name>.py. You may alternatively want to include
# your mod in an existing file, if it coexists naturally with other mods.
#
# In general the public API should be a single class or function. Go with
# whatever makes the most sense for this mod.

import gurobipy as gp
from gurobipy import GRB
import networkx as nx


def solve_network_design(G: nx.DiGraph, commodities):
    """
    A sphinx-compatible docstring

    :param data1: Data structure for first argument
    :type data1: pd.DataFrame
    """
    print(commodities)
    with gp.Env() as env, gp.Model(env=env) as model:

        m = gp.Model()
        x = {
            (i, j, k): m.addVar(name=f"flow-{i},{j},{k}")
            for (i, j) in G.edges()
            for k in commodities
        }
        nodes = G.nodes()
        arcs = G.edges()

        return
        # x = m.addVars(A, C, vtype=GRB.CONTINUOUS, name="flow")
        # y = m.addVars(A, vtype=GRB.BINARY, name="open")

        # m.setObjective(
        #     gp.quicksum(
        #         gp.quicksum(
        #             arc.flow_cost * commodity.demand * x[i, j, k]
        #             for k, commodity in data.commodities.items()
        #         )
        #         + arc.fixed_cost * y[i, j]
        #         for (i, j), arc in G.data.arcs.items()
        #     )
        # )

        # path_constraints = {
        #     (i, k): m.addConstr(
        #         gp.quicksum([x[i, j, k] for j in node.succecessors])
        #         - gp.quicksum([x[j, i, k] for j in node.predecessors])
        #         == (
        #             1
        #             if i == commodity.origin
        #             else (-1 if i == commodity.destination else 0)
        #         ),
        #         name=f"path-{i}-{k}",
        #     )
        #     for i, node in data.nodes.items()
        #     for k, commodity in data.commodities.items()
        # }

        # capacity_constraints = {
        #     (i, j): m.addConstr(
        #         gp.quicksum(
        #             commodity.demand * x[i, j, k]
        #             for k, commodity in data.commodities.items()
        #         )
        #         <= arc.capacity * y[i, j],
        #         name=f"capacity-{i}-{j}",
        #     )
        #     for (i, j), arc in data.arcs.items()
        # }

        # build model
        model.optimize()
        # post-process and return solution
        return None
