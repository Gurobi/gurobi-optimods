import collections

import numpy as np
import gurobipy as gp
from gurobipy import GRB

import networkx as nx


def network_flow(G: nx.DiGraph, Vd=None):
    """
    Network flow problem

    :param G: networkx.DiGraph with edge and node attributes.
    :param source: node of graph
    :param sink: node of graph
    :param demand: float with demand, optional
    :param Vd: subset of nodes to account for in flow constraints. Optional.
    """

    def _get_demand(n):
        if use_demand:
            return G.nodes[n]["demand"]
        else:
            return 0

    with gp.Env() as env, gp.Model(env=env) as model:
        # Check if demand should be used
        use_demand = False
        if len(nx.get_node_attributes(G, "demand")) > 0:
            use_demand = True

        edges, capacities, costs = gp.multidict(
            {(e[0], e[1]): [e[2]["capacity"], e[2]["cost"]] for e in G.edges(data=True)}
        )

        if Vd is None:
            Vd = list(G.nodes())
        x = model.addVars(edges, lb=0, ub=capacities, name="x", obj=costs)
        model.addConstrs(
            (
                gp.quicksum(x[n, j] for j in G.successors(n))
                - gp.quicksum(x[j, n] for j in G.predecessors(n))
                == _get_demand(n)
                for n in Vd
            ),
            name="flow_balance",
        )

        # Use network Simplex
        # TODO: report bug and re-enable later.
        # model.params.NetworkAlg = 1

        model.write("model.lp")

        model.optimize()

        if model.status == GRB.OPTIMAL:
            return [(e, x[e].X) for e in edges if x[e].X > 1e-3], model.ObjVal
        elif model.status in [GRB.INF_OR_UNBD, GRB.INFEASIBLE]:
            raise ValueError("Your problem is infeasible")
        else:
            raise Exception("An error occurred")


def shortest_path(G, source, sink):
    """
    Cost attributes are user-defined

    :param data1: Data structure for first argument
    :type data1: pd.DataFrame
    """
    # Set capacity attribute
    nx.set_edge_attributes(G, 1, "capacity")
    # Set demand attribute
    nx.set_node_attributes(G, 0, "demand")

    G.nodes[source]["demand"] = 1
    G.nodes[sink]["demand"] = -1

    return network_flow(G)


def min_cost_flow(G, source=None, sink=None):
    """
    A sphinx-compatible docstring

    :param data1: Data structure for first argument
    :type data1: pd.DataFrame
    """
    return network_flow(G)


def max_flow(G, source, sink):
    """
    capacity attributes are user-defined so we only define the cost.

    :param data1: Data structure for first argument
    :type data1: pd.DataFrame
    """
    # Set cost attributes to 0 for all edges
    nx.set_edge_attributes(G, 0, "cost")
    # Set demand to 0 for all nodes
    nx.set_node_attributes(G, 0, "demand")
    # Set cost attribute to -1 for all outgoing edges from source
    for j in G.successors(source):
        G.edges[(source, j)]["cost"] = -1
    return network_flow(
        G,
        Vd=[n for n in G.nodes() if n not in [source, sink]],
    )


def min_cut(G, capacities, source, sink):
    """
    A sphinx-compatible docstring

    :param data1: Data structure for first argument
    :type data1: pd.DataFrame
    """
    pass
