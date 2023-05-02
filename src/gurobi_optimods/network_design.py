# Implementation of your new mod. This should be copied to
# src/gurobi_optimods/<mod-name>.py. You may alternatively want to include
# your mod in an existing file, if it coexists naturally with other mods.
#
# In general the public API should be a single class or function. Go with
# whatever makes the most sense for this mod.

import gurobipy as gp
from gurobipy import GRB
import networkx as nx
from typing import Dict


def solve_network_design(G: nx.DiGraph, commodities: Dict):
    """
    A sphinx-compatible docstring

    :param data1: Data structure for first argument
    :type data1: pd.DataFrame
    """
    print(commodities)
    with gp.Env() as env, gp.Model(env=env) as m:

        x = {
            (i, j, k): m.addVar(vtype=GRB.CONTINUOUS, name=f"flow-{i},{j},{k}")
            for (i, j) in G.edges()
            for k in commodities
        }
        y = {
            (i, j): m.addVar(vtype=GRB.BINARY, name=f"arc-{i},{j}")
            for (i, j) in G.edges()
        }

        m.setObjective(
            gp.quicksum(
                gp.quicksum(
                    edge["flow_cost"] * commodity["Demand"] * x[i, j, k]
                    for k, commodity in commodities.items()
                )
                + edge["fixed_cost"] * y[i, j]
                for (i, j, edge) in G.edges(data=True)
            )
        )

        path_constraints = {
            (i, k): m.addConstr(
                gp.quicksum([x[i, j, k] for j in G.successors(i)])
                - gp.quicksum([x[j, i, k] for j in G.predecessors(i)])
                == (
                    1
                    if i == commodity["Origin"]
                    else (-1 if i == commodity["Destination"] else 0)
                ),
                name=f"path-{i}-{k}",
            )
            for i in G.nodes()
            for k, commodity in commodities.items()
        }

        capacity_constraints = {
            (i, j): m.addConstr(
                gp.quicksum(
                    commodity["Demand"] * x[i, j, k]
                    for k, commodity in commodities.items()
                )
                <= arc["capacity"] * y[i, j],
                name=f"capacity-{i}-{j}",
            )
            for (i, j, arc) in G.edges(data=True)
        }
        m.optimize()
        if m.Status == GRB.INFEASIBLE:
            raise ValueError("Unsatisfiable flows")

        G_new = nx.DiGraph()
        for (i, j) in G.edges():
            if y[(i, j)].X > 0.5:
                print(y[(i, j)].X)
                G_new.add_edge(i, j, flow={k: x[i, j, k].X for k in commodities})

        return m.ObjVal, G_new
