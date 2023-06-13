"""
Minimum Cut
-----------
"""

import logging
from typing import Optional, overload

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

from gurobi_optimods.utils import optimod

logger = logging.getLogger(__name__)


@overload
def min_cut(
    graph: sp.spmatrix,
    source: int,
    sink: int,
    silent: bool = False,
    logfile: Optional[str] = None,
) -> sp.spmatrix:
    ...


@overload
def min_cut(
    graph: pd.DataFrame,
    source: int,
    sink: int,
    silent: bool = False,
    logfile: Optional[str] = None,
) -> pd.DataFrame:
    ...


if nx is not None:

    @overload
    def min_cut(
        graph: nx.Graph,
        source: int,
        sink: int,
        silent: bool = False,
        logfile: Optional[str] = None,
    ) -> nx.Graph:
        ...


@optimod()
def min_cut(graph, source: int, sink: int, *, create_env):
    """Solve the minimum cut problem for a given graph.

    :param graph: A graph, specified either as a scipy.sparse adjacency matrix, networkx
        graph, or pandas dataframe
    :type graph: :class:`sp.sparray|nx.Graph|pd.DataFrame`
    :param source: The source node for the path.
    :param sink: The sink (or destination) node for the path.
    :param silent: silent=True suppresses all console output (defaults to False)
    :type silent: bool
    :param logfile: Write all mod output to the given file path (defaults to
        None: no log)
    :type logfile: str
    :return: Cut value of the minimum cut.
    :rtype: :class:`float`
    :return: Partition of size 2 with cut sets.
    :rtype: :class:`tuple(set(), set())`
    :return: Cutset with edges.
    :rtype: :class:`set(tuple(.,.))`
    """
    if isinstance(graph, sp.spmatrix):
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

        if model.Status in [GRB.INFEASIBLE, GRB.INF_OR_UNBD]:
            raise ValueError("Unsatisfiable flows")

        # Construct partition and cutset

        cap_pi = capacity_constrs.gppd.Pi

        cutset = set([a for a in arc_data.index if cap_pi[a] > 1e-3])

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
        return model.ObjVal, (p1, p2), cutset


def _min_cut_scipy(G, source, sink, create_env):
    # Create new matrix with dummy edge (sink, source)
    max_flow = G.tolil()[[0], :].sum()

    data = np.append(G.data, max_flow)

    from_arc = np.append(G.row, sink)
    to_arc = np.append(G.col, source)
    G = sp.coo_matrix((data, (from_arc, to_arc)), dtype=float)

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
        model.update()
        cap = model.addConstr(x <= capacities, name="capacity")
        flow = model.addMConstr(A, x, GRB.EQUAL, demands, name="flow")
        model.optimize()

        if model.Status in [GRB.INFEASIBLE, GRB.INF_OR_UNBD]:
            raise ValueError("Unsatisfiable flows")

        # Construct partition and cutset
        cap_pi = cap.Pi

        p1 = set()
        cutset = set(
            [(i, j) for n, (i, j) in enumerate(zip(G.row, G.col)) if cap_pi[n] > 1e-3]
        )
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
        return model.ObjVal, (p1, p2), cutset


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
        G.remove_edge(sink, source)

        if model.Status in [GRB.INFEASIBLE, GRB.INF_OR_UNBD]:
            raise ValueError("Unsatisfiable flows")

        # Construct partition and cutset
        cutset = set([(i, j) for (i, j) in edges if capacity_constrs[i, j].Pi > 1e-3])

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
        return model.ObjVal, (p1, p2), cutset
