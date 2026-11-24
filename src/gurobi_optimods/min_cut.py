"""
Minimum Cut
-----------
"""

import logging
from dataclasses import dataclass

import gurobipy as gp
import gurobipy_pandas as gppd  # noqa
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
    # Create dummy edge (sink, source) to carry the flow. Its flow variable has
    # no upper bound (hence no bound dual), so it can never be part of the min cut.
    f, t = arc_data.index.names
    arc_data = pd.concat(
        [
            arc_data.assign(cost=0),
            pd.DataFrame(
                [{f: sink, t: source, "capacity": np.inf, "cost": 1}]
            ).set_index([f, t]),
        ]
    )

    with create_env() as env, gp.Model(env=env) as model:
        # Solve max-flow problem
        model.ModelSense = GRB.MAXIMIZE
        arc_df = arc_data.gppd.add_vars(model, obj="cost", name="flow", ub="capacity")
        balance_df = pd.DataFrame(
            {
                "inflow": arc_df["flow"].groupby(t).sum(),
                "outflow": arc_df["flow"].groupby(f).sum(),
            }
        ).gppd.add_constrs(model, "inflow == outflow", name="balance")
        logger.info(
            f"Solving min-cut problem with {len(balance_df)} nodes and "
            f"{len(arc_data)-1} edges"
        )
        model.optimize()

        if model.Status in [GRB.INFEASIBLE, GRB.INF_OR_UNBD]:
            raise ValueError("Unsatisfiable flows")

        # Cutset: arcs whose capacity bound is binding (positive reduced cost).
        # The unbounded dummy (sink, source) arc always has a reduced cost <= 0.
        reduced_costs = arc_df["flow"].gppd.RC
        cutset = set(reduced_costs.index[reduced_costs > 1e-3])
        assert (sink, source) not in cutset

        if len(cutset) == 0 or model.ObjVal < 1e-6:  # No flow, no cut
            return MinCutResult(0.0, (set(), set()), set())

        # Construct partition. The dummy arc is still in arc_data, but the BFS
        # cannot reach the sink, so it is never traversed.
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
        p2 = {n for n in balance_df.index if n not in p1}
        return MinCutResult(model.ObjVal, (p1, p2), cutset)


def _min_cut_scipy(G, source, sink, create_env):
    # Create new matrix with dummy edge (sink, source) to carry the flow. Its
    # flow variable has no upper bound (hence no bound dual), so it can never be
    # part of the min cut.
    data = np.append(G.data, np.inf)

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
        x = model.addMVar(A.shape[1], lb=0, ub=capacities, obj=costs, name="x")
        model.addMConstr(A, x, GRB.EQUAL, demands, name="flow")
        model.optimize()

        if model.Status in [GRB.INFEASIBLE, GRB.INF_OR_UNBD]:
            raise ValueError("Unsatisfiable flows")

        # Cutset: arcs whose capacity bound is binding (positive reduced cost).
        # The unbounded dummy (sink, source) arc (last entry) always has a
        # reduced cost <= 0.
        in_cutset = x.RC > 1e-3
        assert not in_cutset[-1]
        cutset = set(zip(from_arc[in_cutset], to_arc[in_cutset]))

        if len(cutset) == 0 or model.ObjVal < 1e-6:  # No flow, no cut
            return MinCutResult(0.0, (set(), set()), set())

        # Construct partition. The dummy arc is still in G, but the BFS cannot
        # reach the sink, so it is never traversed.
        p1 = set()
        queue = [source]
        G = G.tocsr()
        while len(queue) > 0:
            node = queue.pop()
            p1.add(node)
            row = G[[node]]
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
    # Dummy edge (sink, source) to carry the flow. Its flow variable has no
    # upper bound (hence no bound dual), so it can never be part of the min cut.
    G.add_edge(sink, source, cost=1, capacity=float("inf"))
    with create_env() as env, gp.Model(env=env) as model:
        model.ModelSense = GRB.MAXIMIZE
        edges, capacities, costs = gp.multidict(
            {(i, j): [d["capacity"], d["cost"]] for i, j, d in G.edges(data=True)}
        )
        nodes = list(G.nodes(data=True))
        x = {
            (i, j): model.addVar(
                ub=capacities[i, j], obj=costs[i, j], name=f"flow[{i},{j}]"
            )
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
        model.optimize()
        G.remove_edge(sink, source)

        if model.Status in [GRB.INFEASIBLE, GRB.INF_OR_UNBD]:
            raise ValueError("Unsatisfiable flows")

        # Cutset: arcs whose capacity bound is binding (positive reduced cost).
        # The unbounded dummy (sink, source) arc always has a reduced cost <= 0.
        cutset = set([(i, j) for (i, j) in edges if x[i, j].RC > 1e-3])
        assert (sink, source) not in cutset
        if len(cutset) == 0 or model.ObjVal < 1e-6:  # No flow, no cut
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
