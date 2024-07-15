"""
Maximum Flow
------------
"""

import logging

import numpy as np
import pandas as pd
import scipy.sparse as sp

try:
    import networkx as nx
except ImportError:
    nx = None

from gurobi_optimods.min_cost_flow import (
    min_cost_flow_networkx,
    min_cost_flow_pandas,
    min_cost_flow_scipy,
)

logger = logging.getLogger(__name__)


def max_flow(graph, source: int, sink: int, **kwargs):
    """Solve the maximum flow problem for a given graph.

    Parameters
    ----------
    graph : spmatrix or Graph or DataFrame
        A graph, specified either as a scipy.sparse adjacency matrix, networkx
        graph, or pandas dataframe. These contain the capacity for each edge. In
        the networkx (pandas dataframe) case, we expect the edge attribute
        (column name) ``capacity``. Please see the example in the documentation.
    source : int or str
        The source (or origin) node for the maximum flow.
    sink : int or str
        The sink (or destination) node for the maximum flow.

    Returns
    -------
    flow: float
        Maximum feasible flow through the network.
    subgraph: spmatrix or Graph or DataFrame
        A subgraph of the original graph specifying the flow.
    """
    if sp.issparse(graph):
        return _max_flow_scipy(graph, source, sink, **kwargs)
    elif isinstance(graph, pd.DataFrame):
        return _max_flow_pandas(graph, source, sink, **kwargs)
    elif nx is not None and isinstance(graph, nx.Graph):
        return _max_flow_networkx(graph, source, sink, **kwargs)
    else:
        raise ValueError(f"Unknown graph type: {type(graph)}")


def _remove_dummy_edge(graph, source, sink):
    if sp.issparse(graph):
        graph = graph.tolil()
        graph = graph[:-1, :]
    elif isinstance(graph, pd.Series) or isinstance(graph, pd.DataFrame):
        f, t = ("source", "target")
        graph.drop(graph[(graph[f] == sink) & (graph[t] == source)].index, inplace=True)
    elif nx is not None and isinstance(graph, nx.Graph):
        graph.remove_edge(sink, source)
    return graph


def _max_flow_pandas(arc_data, source, sink, **kwargs):
    f, t = ("source", "target")
    arc_data["cost"] = [0] * len(arc_data)
    # Find maximum flow through (minumum of sum of all outgoing/incoming
    # capacities at the source/sink)
    max_flow = min(
        arc_data[arc_data[f] == source]["capacity"].sum(),
        arc_data[arc_data[t] == sink]["capacity"].sum(),
    )
    # Add dummy edge
    arc_data = pd.concat(
        [
            arc_data,
            pd.DataFrame([{f: sink, t: source, "capacity": max_flow, "cost": -1}]),
        ]
    ).reset_index(drop=True)
    demand_data = pd.DataFrame([{"node": source, "demand": 0}]).set_index("node")
    # Solve
    obj, flow = min_cost_flow_pandas(arc_data, demand_data, **kwargs)
    arc_data = _remove_dummy_edge(arc_data, source, sink)
    flow = _remove_dummy_edge(flow, source, sink)
    return -obj, flow


def _max_flow_scipy(G, source, sink, **kwargs):
    # Create new matrix with dummy edge (sink, source)
    max_flow = G.tolil()[[0], :].sum()
    data = np.append(G.data, max_flow)
    from_arc = np.append(G.row, sink)
    to_arc = np.append(G.col, source)
    G = sp.coo_array((data, (from_arc, to_arc)), dtype=float)

    capacities = sp.coo_array(G)

    costs = np.zeros(G.row.shape, dtype=float)
    costs[-1] = -1
    costs = sp.coo_array((costs, (G.row, G.col)), dtype=float)
    demands = np.zeros(G.shape[1], dtype=float)
    # Solve
    obj, flow = min_cost_flow_scipy(G, capacities, costs, demands, **kwargs)
    G = _remove_dummy_edge(G, source, sink)
    flow = _remove_dummy_edge(flow, source, sink)
    return -obj, sp.coo_array(flow)


def _max_flow_networkx(G, source, sink, **kwargs):
    nx.set_edge_attributes(G, 0, "cost")
    nx.set_node_attributes(G, 0, "demand")
    max_flow = 0
    for j in G.successors(source):
        # G.edges[source, j]["cost"] = -1
        max_flow += G.edges[source, j]["capacity"]
    G.add_edge(sink, source, cost=-1, capacity=max_flow)
    # Set cost attribute to -1 for all outgoing edges from source
    obj, flow_graph = min_cost_flow_networkx(G, **kwargs)
    G = _remove_dummy_edge(G, source, sink)
    if len(flow_graph.edges) > 0:
        flow_graph = _remove_dummy_edge(flow_graph, source, sink)
    return -obj, flow_graph
