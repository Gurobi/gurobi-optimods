"""
Minimum Cut
-----------
"""

import logging
from dataclasses import dataclass

import gurobipy as gp
import gurobipy_pandas as gppd
import numpy as np
import pandas as pd
import scipy.sparse as sp
from gurobipy import GRB

try:
    import networkx as nx
except ImportError:
    nx = None

from gurobi_optimods.max_flow import _remove_dummy_edge
from gurobi_optimods.utils import optimod

logger = logging.getLogger(__name__)


@dataclass
class MinCutResult:
    """
    Solution to a Minimum-cut problem.

    Attributes
    ----------
    cut: float
        Cut value of the minimum cut.
    partition: tuple with sets
        Partition of size 2 with cut sets.
    cutset: set of tuple
        Cutset with edges.
    """

    cut_value: float
    partition: tuple
    cutset: set


@optimod()
def min_cut(graph, source, sink, *, create_env):
    """Solve the minimum cut problem for a given graph.

    Parameters
    ----------
    graph : spmatrix or Graph or DataFrame
        A graph, specified either as a scipy.sparse adjacency matrix, networkx
        graph, or pandas dataframe. These contain the capacity for each edge. In
        the networkx (pandas dataframe) case, we expect the edge attribute
        (column name) ``capacity``. Please see the example in the documentation.
    source : int or str
        The source (or origin) node for the cutset.
    sink : int or str
        The sink (or destination) node for the cutset.

    Returns
    -------
    min_cut_result: MinCutResult
        A dataclass containing the cut value, and set of nodes and edges in the
        minimum cut.
    """
    if sp.issparse(graph):
        return _min_cut_scipy(graph, source, sink, create_env)
    elif isinstance(graph, pd.DataFrame):
        return _min_cut_pandas(graph, source, sink, create_env)
    elif nx is not None and isinstance(graph, nx.Graph):
        return _min_cut_networkx(graph, source, sink, create_env)
    else:
        raise ValueError(f"Unknown graph type: {type(graph)}")


def _min_cut_pandas(arc_data, source, sink, create_env):
    f, t = arc_data.index.names
    arc_data["cost"] = [0] * len(arc_data)
    # Create dummy edge to find maximum flow through (minimum of sum of all
    # outgoing/incoming capacities at the source/sink)
    max_flow = min(
        arc_data.loc[([source], slice(None)),]["capacity"].sum(),
        arc_data.loc[(slice(None), [sink]),]["capacity"].sum(),
    )
    arc_data = pd.concat(
        [
            arc_data,
            pd.DataFrame(
                [{f: sink, t: source, "capacity": max_flow, "cost": 1}]
            ).set_index([f, t]),
        ]
    )

    with create_env() as env, gp.Model(env=env) as model:
        # Solve max-flow problem
        model.ModelSense = GRB.MAXIMIZE
        arc_df = arc_data.gppd.add_vars(model, obj="cost", name="flow")
        balance_df = pd.DataFrame(
            {
                "inflow": arc_df["flow"].groupby(t).sum(),
                "outflow": arc_df["flow"].groupby(f).sum(),
            }
        ).gppd.add_constrs(model, "inflow == outflow", name="balance")
        capacity_constrs = gppd.add_constrs(
            model,
            arc_df["flow"],
            GRB.LESS_EQUAL,
            arc_df["capacity"],
            name="capacity",
        )
        logger.info(
            f"Solving min-cut problem with {len(balance_df)} nodes and "
            f"{len(arc_data)-1} edges"
        )
        model.optimize()
        _remove_dummy_edge(arc_data, source, sink)

        if model.Status in [GRB.INFEASIBLE, GRB.INF_OR_UNBD]:
            raise ValueError("Unsatisfiable flows")

        # Construct partition and cutset
        cap_pi = capacity_constrs.gppd.Pi

        # Find edges in the cutset (excluding the dummy (sink, source) edge.
        cutset = set([a for a in arc_data.index if cap_pi[a] > 1e-3 if a[0] != sink])

        if len(cutset) == 0:  # No arc in the cutset
            return MinCutResult(0.0, (set(), set()), set())

        p1 = set()
        queue = [source]
        while len(queue) > 0:
            node = queue.pop()
            p1.add(node)
            # Add successors of `node` that are not in the cutset
            queue.extend(
                [
                    a[1]
                    for a in arc_data.loc[([node], slice(None)),].index
                    if a not in cutset and a[1] not in p1 and a[1] not in queue
                ]
            )
        p2 = set([n for n in balance_df.index if n not in p1])
        return MinCutResult(model.ObjVal, (p1, p2), cutset)


def _min_cut_scipy(G, source, sink, create_env):
    # Create new matrix with dummy edge (sink, source)
    original_shape = G.shape
    max_flow = G.tolil()[[0], :].sum()

    data = np.append(G.data, max_flow)

    from_arc = np.append(G.row, sink)
    to_arc = np.append(G.col, source)
    G = sp.coo_array((data, (from_arc, to_arc)), dtype=float)

    capacities = data

    costs = np.zeros(G.row.shape, dtype=float)
    costs[-1] = 1
    demands = np.zeros(G.shape[1], dtype=float)

    G = G.tocoo()

    # Create incidence matrix from edge lists.
    indices = np.column_stack((from_arc, to_arc)).reshape(-1, order="C")
    indptr = np.arange(0, 2 * from_arc.shape[0] + 2, 2)
    ones = np.ones(from_arc.shape)
    data = np.column_stack((ones * -1.0, ones)).reshape(-1, order="C")

    A = sp.csc_array((data, indices, indptr))

    logger.info(
        f"Solving min-cut problem with {A.shape[0]} nodes and " f"{A.shape[1]-1} edges"
    )

    with create_env() as env, gp.Model(env=env) as model:
        # Solve max-flow problem
        model.ModelSense = GRB.MAXIMIZE
        x = model.addMVar(A.shape[1], lb=0, obj=costs, name="x")
        cap = model.addConstr(x <= capacities, name="capacity")
        model.addMConstr(A, x, GRB.EQUAL, demands, name="flow")
        model.optimize()
        _remove_dummy_edge(G, source, sink)

        if model.Status in [GRB.INFEASIBLE, GRB.INF_OR_UNBD]:
            raise ValueError("Unsatisfiable flows")

        # Construct partition and cutset
        cap_pi = cap.Pi

        p1 = set()
        cutset = set(
            [
                (i, j)
                for n, (i, j) in enumerate(zip(G.row, G.col))
                if cap_pi[n] > 1e-3 and i != sink
            ]
        )
        if len(cutset) == 0:  # No arc in the cutset
            return MinCutResult(0.0, (set(), set()), set())

        queue = [source]
        while len(queue) > 0:
            node = queue.pop()
            p1.add(node)
            row = G.getrow(node)
            # Add successors of `node` that are not in the cutset
            queue.extend(
                [
                    j
                    for j in row.tocoo().col
                    if (node, j) not in cutset and j not in p1 and j not in queue
                ]
            )
        p2 = set([n for n in range(G.shape[1]) if n not in p1])
        return MinCutResult(model.ObjVal, (p1, p2), cutset)


def _min_cut_networkx(G, source, sink, create_env):
    logger.info(
        f"Solving min-cut problem with {len(G.nodes)} nodes and "
        f"{len(G.edges)} edges"
    )
    nx.set_edge_attributes(G, 0, "cost")
    max_flow = 0
    for j in G.successors(source):
        max_flow += G.edges[source, j]["capacity"]
    G.add_edge(sink, source, cost=1, capacity=max_flow)
    with create_env() as env, gp.Model(env=env) as model:
        model.ModelSense = GRB.MAXIMIZE
        edges, capacities, costs = gp.multidict(
            {(i, j): [d["capacity"], d["cost"]] for i, j, d in G.edges(data=True)}
        )
        nodes = list(G.nodes(data=True))
        x = {
            (i, j): model.addVar(obj=costs[i, j], name=f"flow[{i},{j}]")
            for i, j in edges
        }

        flow_constrs = {
            n: model.addConstr(
                (
                    gp.quicksum(x[j, n] for j in G.predecessors(n))
                    == gp.quicksum(x[n, j] for j in G.successors(n))
                ),
                name=f"flow_balance[{n}]",
            )
            for n, data in nodes
        }
        model.update()
        capacity_constrs = {
            (i, j): model.addConstr(
                x[i, j] <= capacities[i, j], name=f"capacity[{i},{j}]"
            )
            for i, j in edges
        }

        model.optimize()
        _remove_dummy_edge(G, source, sink)

        if model.Status in [GRB.INFEASIBLE, GRB.INF_OR_UNBD]:
            raise ValueError("Unsatisfiable flows")

        # Construct partition and cutset
        cutset = set(
            [
                (i, j)
                for (i, j) in edges
                if capacity_constrs[i, j].Pi > 1e-3 and i != sink
            ]
        )
        if len(cutset) == 0:
            return MinCutResult(0.0, (set(), set()), set())

        p1 = set()
        queue = [source]
        while len(queue) > 0:
            node = queue.pop()
            p1.add(node)
            # Add successors of `node` that are not in the cutset
            queue.extend(
                [
                    j
                    for j in G.successors(node)
                    if (node, j) not in cutset and j not in p1 and j not in queue
                ]
            )
        p2 = set([n for n in G.nodes if n not in p1])
        return MinCutResult(model.ObjVal, (p1, p2), cutset)
