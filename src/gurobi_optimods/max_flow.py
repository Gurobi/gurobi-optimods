import logging
from typing import Optional, overload

import numpy as np
import pandas as pd
import scipy.sparse as sp

try:
    import networkx as nx
except ImportError:
    nx = None

from gurobi_optimods.utils import optimod
from gurobi_optimods.min_cost_flow import (
    min_cost_flow,
    min_cost_flow_scipy,
    min_cost_flow_networkx,
)

logger = logging.getLogger(__name__)


@overload
def max_flow(
    graph: sp.spmatrix,
    source: int,
    sink: int,
    silent: bool = False,
    logfile: Optional[str] = None,
) -> sp.spmatrix:
    ...


@overload
def max_flow(
    graph: pd.DataFrame,
    source: int,
    sink: int,
    silent: bool = False,
    logfile: Optional[str] = None,
) -> pd.DataFrame:
    ...


if nx is not None:

    @overload
    def max_flow(
        graph: nx.Graph,
        source: int,
        sink: int,
        silent: bool = False,
        logfile: Optional[str] = None,
    ) -> nx.Graph:
        ...


def max_flow(graph, source: int, sink: int, **kwargs):
    """Solve the maximum flow problem for a given graph.

    :param graph: A graph, specified either as a scipy.sparse adjacency matrix, networkx
        graph, or pandas dataframe
    :type graph: :class:`sp.sparray|nx.Graph|pd.DataFrame`
    :param source: The source node for the path.
    :type source: :class:`int`
    :param sink: The sink (or destination) node for the path.
    :type sink: :class:`int`
    :param silent: silent=True suppresses all console output (defaults to False)
    :type silent: bool
    :param logfile: Write all mod output to the given file path (defaults to None: no log)
    :type logfile: str
    :return: Maximum feasible flow through the network.
    :rtype: :class:`float`
    :return: A subgraph of the original graph specifying the flow.
    :rtype: :class:`sp.sparray|nx.Graph|pd.DataFrame`
    """
    if isinstance(graph, sp.spmatrix):
        return _max_flow_scipy(graph, source, sink, **kwargs)
    elif isinstance(graph, pd.DataFrame):
        return _max_flow_pandas(graph, source, sink, **kwargs)
    elif nx is not None and isinstance(graph, nx.Graph):
        return _max_flow_networkx(graph, source, sink, **kwargs)
    else:
        raise ValueError(f"Unknown graph type: {type(graph)}")


def _max_flow_pandas(arc_data, source, sink, **kwargs):
    f, t = arc_data.index.names
    arc_data["cost"] = [0] * len(arc_data)
    # Find maximum flow through (minumum of sum of all outgoing/incoming
    # capacities at the source/sink)
    max_flow = min(
        arc_data.loc[
            ([source], slice(None)),
        ]["capacity"].sum(),
        arc_data.loc[
            (slice(None), [sink]),
        ]["capacity"].sum(),
    )
    # Add dummy edge
    arc_data = pd.concat(
        [
            arc_data,
            pd.DataFrame(
                [{f: sink, t: source, "capacity": max_flow, "cost": -1}]
            ).set_index([f, t]),
        ]
    )
    demand_data = pd.DataFrame([{"node": source, "demand": 0}]).set_index("node")
    # Solve
    obj, flow = min_cost_flow(arc_data, demand_data, **kwargs)
    # Drop dummy edge from solution
    flow.drop((sink, source), inplace=True)
    return -obj, flow


def _max_flow_scipy(G, source, sink, **kwargs):
    # Create new matrix with dummy edge (sink, source)
    original_shape = G.shape
    max_flow = G.tolil()[[0], :].sum()
    data = np.append(G.data, max_flow)
    from_arc = np.append(G.row, sink)
    to_arc = np.append(G.col, source)
    G = sp.coo_matrix((data, (from_arc, to_arc)), dtype=float)

    capacities = sp.coo_matrix(G)

    costs = np.zeros(G.row.shape, dtype=float)
    costs[-1] = -1
    costs = sp.coo_matrix((costs, (G.row, G.col)), dtype=float)
    demands = np.zeros(G.shape[1], dtype=float)
    # Solve
    obj, flow = min_cost_flow_scipy(G, capacities, costs, demands, **kwargs)
    # Remove dummy edge
    flow = flow.tolil()
    # Remove row
    if flow.shape != original_shape:
        flow = flow[:-1, :]
    return -obj, sp.coo_matrix(flow)


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
    flow_graph.remove_edge(sink, source)
    G.remove_edge(sink, source)
    return -obj, flow_graph
