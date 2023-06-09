"""
Bipartite Matching
------------------
"""


import logging

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


@optimod()
def maximum_bipartite_matching(graph, nodes1, nodes2, *, create_env):
    """Solve a maximum cardinality bipartite matching problem on the
    given graph.

    Parameters
    ----------
    graph : spmatrix or Graph or DataFrame
        A graph, specified either as a scipy.sparse adjacency matrix, networkx
        graph, or pandas dataframe.
    nodes1 : ndarray or str
        Nodes in the first bipartite set. If ``graph`` is a scipy sparse matrix,
        ``nodes1`` must be a numpy array. If ``graph`` is a pandas dataframe,
        ``nodes1`` must be a column name. If ``graph`` is a networkx graph,
        ``nodes1`` must be a list.
    nodes2 : ndarray or str
        Nodes in the first bipartite set. If ``graph`` is a scipy sparse matrix,
        ``nodes2`` must be a numpy array. If ``graph`` is a pandas dataframe,
        ``nodes2`` must be a column name. If ``graph`` is a networkx graph,
        ``nodes2`` must be a list.

    Returns
    -------
    spmatrix or Graph or DataFrame
        A subgraph of the original ``graph`` (with the same data type) specifying
        the maximum matching
    """
    if isinstance(graph, sp.spmatrix):
        return _maximum_bipartite_matching_scipy(graph, nodes1, nodes2, create_env)
    elif isinstance(graph, pd.DataFrame):
        return _maximum_bipartite_matching_pandas(graph, nodes1, nodes2, create_env)
    elif nx is not None and isinstance(graph, nx.Graph):
        return _maximum_bipartite_matching_networkx(graph, nodes1, nodes2, create_env)
    else:
        raise ValueError(f"Unknown graph type: {type(graph)}")


def _maximum_bipartite_matching_pandas(frame, n1_column, n2_column, create_env):
    """This implementation uses gurobipy-pandas, which suits the input data
    already in a pandas dataframe."""

    with create_env() as env, gp.Model(env=env) as model:
        # Directed flow variables between bipartite sets
        df = frame.set_index([n1_column, n2_column]).gppd.add_vars(
            model, ub=1, name="flow"
        )

        # Inflow variables and flow balance on n1
        n1_outflows = df["flow"].groupby(n1_column).sum()
        n1_inflows = gppd.add_vars(model, n1_outflows.index, ub=1, name="src_flow")
        gppd.add_constrs(model, n1_inflows, GRB.EQUAL, n1_outflows, name="n1_balance")

        # Outflow variables and flow balance on n2
        n2_inflows = df["flow"].groupby(n2_column).sum()
        n2_outflows = gppd.add_vars(model, n2_inflows.index, ub=1, name="sink_flow")
        gppd.add_constrs(model, n2_inflows, GRB.EQUAL, n2_outflows, name="n2_balance")

        # sink-source flow variable and flow balances
        sink_source_flow = model.addVar(name="sink_source_flow")
        model.addConstr(sink_source_flow == n1_inflows.sum(), name="src_balance")
        model.addConstr(n2_outflows.sum() == sink_source_flow, name="sink_balance")

        # max flow
        model.setObjective(sink_source_flow, sense=GRB.MAXIMIZE)

        # solve and extract solution
        model.optimize()
        return df.loc[df["flow"].gppd.X.gt(0.1)].reset_index().drop(columns=["flow"])


def _maximum_bipartite_matching_networkx(graph, nodes1, nodes2, create_env):
    """This implementation uses gurobipy's term-based API, which suits the
    iterator-based API for reading data from networkx graphs."""

    logger.info(
        f"Solving maximum matching n1={len(nodes1)} "
        f"n2={len(nodes2)} |E|={graph.number_of_edges()}"
    )

    nodes1 = set(nodes1)
    nodes2 = set(nodes2)

    # Treat all matching problems as undirected
    graph = graph.to_undirected()

    with create_env() as env, gp.Model(env=env) as model:
        # Add variables for each layer of edges in the max flow graph
        source_layer = {i: model.addVar(name=f"flow[source,{i}]", ub=1) for i in nodes1}
        graph_layer = {
            (i, j): model.addVar(name=f"flow[{i},{j}]", ub=1)
            for i, j in graph.edges
            if i in nodes1 and j in nodes2
        }
        sink_layer = {j: model.addVar(name=f"flow[{j},sink]", ub=1) for j in nodes2}
        sink_source = model.addVar(name="flow[sink,source]")

        # At the source node, sink -> source flow balances source -> nodes1 flow
        model.addConstr(
            gp.quicksum(source_layer.values()) == sink_source, name="source_balance"
        )

        # In the nodes1 layer, flows from source balance flows to nodes2
        for i in nodes1:
            model.addConstr(
                source_layer[i] == gp.quicksum(graph_layer[i, j] for j in graph[i]),
                name=f"n1_balance[{i}]",
            )

        # In the nodes1 layer, flows from nodes2 balance flows to sink
        for j in nodes2:
            model.addConstr(
                sink_layer[j] == gp.quicksum(graph_layer[i, j] for i in graph[j]),
                name=f"n2_balance[{j}]",
            )

        # At the sink node, sink -> source flow balances nodes2 -> sink flow
        model.addConstr(
            gp.quicksum(sink_layer.values()) == sink_source, name="sink_balance"
        )

        # Maximize flow through the uncapacitated sink->source edge
        model.setObjective(sink_source, sense=GRB.MAXIMIZE)
        model.optimize()

        # Create a new Graph with selected edges in the matching
        matching = nx.Graph()
        matching.add_nodes_from(graph.nodes)
        matching.add_edges_from(edge for edge, x in graph_layer.items() if x.X > 0.1)

        logger.info(f"Max bipartite matching |E|={matching.number_of_edges()}")

        return matching


def _maximum_bipartite_matching_scipy(adjacency, nodes1, nodes2, create_env):
    """This implementation uses the gurobipy matrix-friendly API, which suits
    the input data in numpy/scipy data structures."""

    logger.info(
        f"Solving maximum matching n1={nodes1.shape[0]} "
        f"n2={nodes2.shape[0]} |E|={adjacency.data.shape[0]}"
    )

    # Add a source and sink node for max flow formulation
    # Assume G is symmetric (or upper triangular)
    G = sp.triu(adjacency.tocoo())
    G_nodes = adjacency.shape[0]
    source, sink = G_nodes, G_nodes + 1

    # Build network:
    #   source -> nodes1 (complete)
    #   nodes1 -> nodes2 (adjacency)
    #   nodes2 -> sink (complete)
    #   sink -> source
    from_arc = np.concatenate([np.repeat(source, nodes1.shape), G.row, nodes2, [sink]])
    to_arc = np.concatenate([nodes1, G.col, np.repeat(sink, nodes2.shape), [source]])
    capacity = np.ones(from_arc.shape, dtype=float)
    capacity[-1] = GRB.INFINITY
    balance = np.zeros(G_nodes + 2)
    logger.info(
        "Maximum matching formulated as min-cost flow with "
        f"{balance.shape[0]} nodes and {from_arc.shape[0]} arcs"
    )

    # Solve min-cost flow problem (in this case, it's actually max-profit)
    with create_env() as env, gp.Model(env=env) as model:
        # Create incidence matrix from edge lists.
        indices = np.column_stack((from_arc, to_arc)).reshape(-1, order="C")
        indptr = np.arange(0, 2 * from_arc.shape[0] + 2, 2)
        ones = np.ones(from_arc.shape)
        data = np.column_stack((ones * -1.0, ones)).reshape(-1, order="C")
        A = sp.csc_array((data, indices, indptr))

        # Solve model with gurobi, return cost and flows
        x = model.addMVar(A.shape[1], lb=0, ub=capacity)
        model.setObjective(x[-1], sense=GRB.MAXIMIZE)
        model.addMConstr(A, x, GRB.EQUAL, balance)
        model.optimize()
        if model.Status == GRB.INFEASIBLE:
            raise ValueError("Unsatisfiable flows")
        flows = x.X

    # Choose the arcs corresponding to the original graph with non-zero
    # flow. Note that the last var is the sink->source connection (drop it).
    select = (flows > 0.5) & (from_arc != source) & (to_arc != sink)
    from_arc_result = from_arc[select][:-1]
    to_arc_result = to_arc[select][:-1]

    logger.info(f"Done: max bipartite matching has {from_arc_result.shape[0]} edges")

    # Return undirected, unweighted adjacency matrix
    arg = (np.ones(from_arc_result.shape), (from_arc_result, to_arc_result))
    matching = sp.coo_matrix(arg, dtype=float, shape=G.shape)
    return matching + matching.T
