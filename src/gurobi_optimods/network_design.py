import gurobipy as gp
from gurobipy import GRB
import networkx as nx
from typing import Dict
from gurobi_optimods.utils import optimod
import matplotlib.pyplot as plt
from gurobi_optimods.datasets import load_network_design


@optimod()
def solve_network_design(G: nx.DiGraph, commodities: Dict, *, create_env):
    """Solve a fixed-cost network design problem on a given graph for a given set of commodities

    :param G: A graph specified as a networkx graph
    :type G: :class:`nx.DiGraph`
    :param commodities: A dictionary where the keys correpond to different commodity labels and the values are triples
        containing the origin node, destination node, and the quantity.
    :type G: :class:`Dict`
    :param silent: silent=True suppresses all console output (defaults to False)
    :type silent: bool
    :param logfile: Write all mod output to the given file path (defaults to None: no log)
    :type logfile: str
    :return: A subgraph of the original graph specifying the maximum matching
    :rtype: :class:`nx.Graph`
    """

    if isinstance(G, nx.Graph):
        return _network_design_networkx(G, commodities, create_env)
    else:
        raise ValueError(f"Unknown graph type: {type(G)}")


def _network_design_networkx(G, commodities, create_env):

    _ensure_origins_desinations_are_nodes(G.nodes(), commodities)
    with create_env() as env, gp.Model(env=env) as m:

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

        m._path_constraints = {
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

        m._capacity_constraints = {
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

        # Create a new graph with selected edges in the match
        G_new = nx.DiGraph()
        for (i, j) in G.edges():
            if y[(i, j)].X > 0.5:
                print(y[(i, j)].X)
                G_new.add_edge(i, j, flow={k: x[i, j, k].X for k in commodities})

        return m.ObjVal, G_new


def _ensure_origins_desinations_are_nodes(nodes, commodities):

    for k, commodity in commodities.items():
        assert (
            commodity["Origin"] in nodes
        ), f"Origin {commodity['Origin']} of commodity {k} not a node"
        assert (
            commodity["Destination"] in nodes
        ), f"Destination {commodity['Destination']} of commodity {k} not a node"


def _draw_network_design_with_solution(original_graph, solution_graph):
    print("hi")

    nx.draw(original_graph, with_labels=True)
    plt.draw()  # pyplot draw()
    plt.show()
