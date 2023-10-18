"""
Minimum Cost Flow
-----------------
"""

import logging

import gurobipy as gp
import gurobipy_pandas as gppd  # noqa: F401
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


@optimod()
def min_cost_flow_pandas(
    arc_data: pd.DataFrame, demand_data: pd.DataFrame, *, create_env
):
    """Solve the minimum cost flow problem for a given graph.

    The inputs adhere to the following structure::

        arc_data = pd.DataFrame(
            [
                {"from": 0, "to": 1, "capacity": 16, "cost": 0}, {"from": 1,
                "to": 2, "capacity": 10, "cost": 0},
            ]
        ).set_index(["from", "to"])

        demand_data = pd.DataFrame(
            [{"node": 0, "demand": -1}, {"node": 2, "demand": 1}]
        ).set_index("node")

    Parameters
    ----------
    arc_data : DataFrame
        DataFrame with graph and respective attributes. These must include
        ``"from"``, ``"to"`` nodes used as index, ``"capacity"``, and
        ``"cost"``.
    demand_data : DataFrame
        DataFrame with node demand information. These must
        include indexed by `"node"`, and include the `"demand"`. This value can
        be positive (requesting flow) or negative (supplying flow).

    Returns
    -------
    tuple
        Cost of the minimum cost flow (float), dictionary indexed by edges with
        non-zero flow in the solution (Series)
    """
    with create_env() as env, gp.Model(env=env) as model:
        model.ModelSense = GRB.MINIMIZE

        arc_df = arc_data.gppd.add_vars(model, ub="capacity", obj="cost", name="flow")

        source_label, target_label = arc_data.index.names
        balance_df = (
            pd.DataFrame(
                {
                    "inflow": arc_df["flow"].groupby(target_label).sum(),
                    "outflow": arc_df["flow"].groupby(source_label).sum(),
                    "demand": demand_data["demand"],
                }
            )
            .fillna(0)  # zero fill (some nodes have no in, out, or demand)
            .gppd.add_constrs(model, "inflow - outflow == demand", name="balance")
        )
        logger.info(
            f"Solving min-cost flow with {len(balance_df)} nodes and "
            f"{len(arc_data)} edges"
        )
        model.optimize()

        if model.Status in [GRB.INFEASIBLE, GRB.INF_OR_UNBD]:
            raise ValueError("Unsatisfiable flows")

        return model.ObjVal, arc_df["flow"].gppd.X


@optimod()
def min_cost_flow_scipy(
    G: sp.spmatrix,
    capacities: sp.spmatrix,
    costs: sp.spmatrix,
    demands: np.ndarray,
    *,
    create_env,
):
    """Solve the minimum cost flow problem for a given graph.

    Parameters
    ----------
    G : spmatrix
        Adjacency matrix of the graph.
    capacities : spmatrix
        Matrix containing capacities for each edge.
    costs : spmatrix
        Matrix containing costs for each edge.
    demands : ndarray
        Array containing the demand for each node.

    Returns
    -------
    tuple
        Cost of the minimum cost flow (float), dictionary indexed by edges with
        non-zero flow in the solution (spmatrix)
    """
    G = G.tocoo()

    edge_source = G.row
    edge_target = G.col

    capacities = capacities.tocoo()
    capacities = capacities.data

    costs = costs.tocoo()
    costs = costs.data

    # Create incidence matrix from edge lists.
    indices = np.column_stack((edge_source, edge_target)).reshape(-1, order="C")
    indptr = np.arange(0, 2 * edge_source.shape[0] + 2, 2)
    ones = np.ones(edge_source.shape)
    data = np.column_stack((ones * -1.0, ones)).reshape(-1, order="C")

    A = sp.csc_array((data, indices, indptr))

    logger.info("Solving min-cost flow with {0} nodes and {1} edges".format(*A.shape))

    # Solve model with gurobi, return cost and flows
    with create_env() as env, gp.Model(env=env) as model:
        x = model.addMVar(A.shape[1], lb=0, obj=costs, name="x")
        model.addConstr(x <= capacities, name="capacity")
        model.addMConstr(A, x, GRB.EQUAL, demands, name="flow")
        model.optimize()
        if model.Status in [GRB.INFEASIBLE, GRB.INF_OR_UNBD]:
            raise ValueError("Unsatisfiable flows")
        # Filter + create scipy output matrix
        select = x.X > 0.5
        edge_source_result = edge_source[select]
        edge_target_result = edge_target[select]
        arg = (x.X[select], (edge_source_result, edge_target_result))
        return model.ObjVal, sp.coo_array(arg, dtype=float, shape=G.shape)


@optimod()
def min_cost_flow_networkx(G, *, create_env):
    """Solve the minimum cost flow problem for a given graph.

    Parameters
    ----------
    G : DiGraph
        Graph with edge attributes ``capacity`` and ``cost``, as well as node
        attributes ``demand``.

    Returns
    -------
    tuple
        Cost of the minimum cost flow (float), a subgraph of the original graph
        specifying the flow
    """
    logger.info(
        f"Solving min-cost flow with {len(G.nodes)} nodes and {len(G.edges)} edges"
    )

    with create_env() as env, gp.Model(env=env) as model:
        edges, capacities, costs = gp.multidict(
            {(i, j): [d["capacity"], d["cost"]] for i, j, d in G.edges(data=True)}
        )
        nodes = list(G.nodes(data=True))
        x = {
            (i, j): model.addVar(
                name=f"flow[{i},{j}]",
                ub=capacities[i, j],
                obj=costs[i, j],
            )
            for i, j in edges
        }

        flow_constrs = {}
        for n, data in nodes:
            flow_constrs[n] = model.addConstr(
                (
                    gp.quicksum(x[j, n] for j in G.predecessors(n))
                    - gp.quicksum(x[n, j] for j in G.successors(n))
                    == data["demand"]
                ),
                name=f"flow_balance[{n}]",
            )

        model.optimize()

        if model.Status in [GRB.INFEASIBLE, GRB.INF_OR_UNBD]:
            raise ValueError("Unsatisfiable flows")

        # Create a new Graph with selected edges in the matching
        resulting_flow = nx.DiGraph()
        resulting_flow.add_nodes_from(nodes)
        resulting_flow.add_edges_from(
            [(edge[0], edge[1], {"flow": v.X}) for edge, v in x.items() if v.X > 0.1]
        )

        return model.ObjVal, resulting_flow
