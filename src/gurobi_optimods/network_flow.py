import numpy as np
import gurobipy as gp
import networkx as nx
from gurobipy import GRB


def network_flow(G: nx.DiGraph, Vd=None, get_partition=False):
    """
    Network flow problem.

    :param G: networkx.DiGraph with edge attributes (capacity and cost) and node attributes (demand).
    :param Vd: subset of nodes to account for in flow constraints. Default is all nodes. Optional.
    :param get_partition: bool whether the min cut partition should be returned. Optional.
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
            {(i, j): [d["capacity"], d["cost"]] for i, j, d in G.edges(data=True)}
        )

        if Vd is None:
            Vd = list(G.nodes())
        x = model.addVars(edges, lb=0, name="x", obj=costs)
        flow_constrs = model.addConstrs(
            (
                gp.quicksum(x[n, j] for j in G.successors(n))
                - gp.quicksum(x[j, n] for j in G.predecessors(n))
                == _get_demand(n)
                for n in Vd
            ),
            name="flow_balance",
        )
        cap_constrs = model.addConstrs(
            (x[i, j] <= capacities[i, j] for (i, j) in G.edges()), name="capacity"
        )

        # Use network Simplex
        # model.params.NetworkAlg = 1
        model.optimize()

        if model.status == GRB.OPTIMAL:
            if get_partition:

                p1 = set({i for (i, j) in edges if cap_constrs[i, j].Pi < -1e-3})
                p2 = set()
                for n in Vd:
                    if flow_constrs[n].Pi < 1e-3:
                        p2.add(n)
                    elif flow_constrs[n].Pi > 1e-3:
                        p1.add(n)
                return (
                    {e: x[e].X for e in edges if x[e].X > 1e-3},
                    model.ObjVal,
                    (p1, p2),
                )
            return {e: x[e].X for e in edges if x[e].X > 1e-3}, model.ObjVal
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


def max_flow(G, source, sink, **kwargs):
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
    if "get_partition" in kwargs and kwargs["get_partition"]:
        sol, flow, partition = network_flow(
            G, Vd=[n for n in G.nodes() if n not in [source, sink]], **kwargs
        )
        return sol, -flow, partition
    sol, flow = network_flow(G, Vd=[n for n in G.nodes() if n not in [source, sink]])
    return sol, -flow


def min_cut(G, source, sink):
    """
    A sphinx-compatible docstring

    :param data1: Data structure for first argument
    :type data1: pd.DataFrame
    """
    r1, r2 = nx.minimum_cut(G, source, sink)
    result, flow, partition = max_flow(G, source, sink, get_partition=True)
    partition[1].add(sink)
    return partition, flow


def read_dimacs_graph(file_path):
    """
    Parse .col file and return graph object
    """
    G = nx.DiGraph()
    with open(file_path, "r") as file:
        for line in file:
            if line.startswith("c"):  # graph description
                print(*line.split()[1:])
            # first line: p name num_of_vertices num_of_edges
            elif line.startswith("p"):
                p, name, vertices_num, edges_num = line.split()
            elif line.startswith("n"):
                _, n, demand = line.split()
                G.add_node(n, demand=float(demand))
            elif line.startswith("e"):
                _, v1, v2 = line.split()
                G.add_edge(v1, v2)
            elif line.startswith("a"):
                _, v1, v2, _, cap, cost = line.split()
                G.add_edge(v1, v2, capacity=float(cap), cost=float(cost))
            else:
                continue
        for n, d in G.nodes(data=True):
            if "demand" not in d:
                G.nodes[n]["demand"] = 0
        return G
