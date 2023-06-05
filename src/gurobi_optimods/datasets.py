"""
Module for loading datasets for use in optimods examples, in the same vein
as sklearn.datasets.
"""

import pathlib

import numpy as np
import pandas as pd
import scipy.sparse as sp

try:
    import networkx as nx
except ImportError:
    nx = None

DATA_FILE_DIR = pathlib.Path(__file__).parent / "data"


class AttrDict(dict):
    """Even simpler version of sklearn's Bunch"""

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(key)


def load_workforce():
    return AttrDict(
        availability=pd.read_csv(
            DATA_FILE_DIR / "workforce/preferences.csv", parse_dates=["Shift"]
        ),
        shift_requirements=pd.read_csv(
            DATA_FILE_DIR / "workforce/shift_requirements.csv", parse_dates=["Shift"]
        ),
        worker_limits=pd.read_csv(
            DATA_FILE_DIR / "workforce/worker_limits.csv",
        ),
    )


def load_graph(drop_pos=True, capacity=True, cost=True, demand=True):
    edge_data = pd.read_csv(DATA_FILE_DIR / "graphs/edge_data1.csv").set_index(
        ["source", "target"]
    )
    node_data = pd.read_csv(DATA_FILE_DIR / "graphs/node_data1.csv", index_col=0)
    if drop_pos:
        node_data.drop(columns=["posx", "posy"], inplace=True)
    if not capacity:
        edge_data.drop(columns=["capacity"], inplace=True)
    if not cost:
        edge_data.drop(columns=["cost"], inplace=True)
    if not demand:
        node_data.drop(columns=["demand"], inplace=True)
    return edge_data, node_data


def load_graph_networkx(capacity=True, cost=True, demand=True):
    edge_data, node_data = load_graph(capacity=capacity, cost=cost, demand=demand)
    return _convert_pandas_to_digraph(
        edge_data, node_data, capacity=capacity, cost=cost, demand=demand
    )


def load_graph_scipy(capacity=True, cost=True, demand=True):
    edge_data, node_data = load_graph(capacity=capacity, cost=cost, demand=demand)
    return _convert_pandas_to_scipy(
        edge_data, node_data, capacity=capacity, cost=cost, demand=demand
    )


def load_graph2():
    return (
        pd.read_csv(DATA_FILE_DIR / "graphs/edge_data2.csv").set_index(
            ["source", "target"]
        ),
        pd.read_csv(DATA_FILE_DIR / "graphs/node_data2.csv", index_col=0),
    )


def load_graph2_networkx():
    edge_data, node_data = load_graph2()
    return _convert_pandas_to_digraph(edge_data, node_data)


def load_graph2_scipy():
    edge_data, node_data = load_graph2()
    return _convert_pandas_to_scipy(edge_data, node_data)


def load_portfolio():
    fn = DATA_FILE_DIR / "portfolio/portfolio.csv"
    return pd.read_csv(fn, index_col=0)


def _convert_pandas_to_digraph(
    edge_data, node_data, capacity=True, cost=True, demand=True
):
    """
    Convert from a pandas DataFrame to a networkx.DiGraph with the appropriate
    attributes. For edges: `capacity`, and `cost`. For nodes: `demand`.
    """
    G = nx.from_pandas_edgelist(
        edge_data.reset_index(), create_using=nx.DiGraph(), edge_attr=True
    )
    if demand:
        for i, d in node_data.iterrows():
            G.add_node(i, demand=d.demand)
    return G


def _convert_pandas_to_scipy(
    edge_data, node_data, capacity=True, cost=True, demand=True
):
    """
    Convert from a pandas DataFrame to several scipy.sparse.coo_matrix contain
    the graph structure, the capacity and cost values per edge, and the demand
    values per node.
    """
    coords = edge_data.index.to_numpy()

    a0 = np.array([c[0] for c in coords])
    a1 = np.array([c[1] for c in coords])

    data = np.ones(len(coords), dtype=np.int64)
    G = sp.coo_matrix((data, (a0, a1)))

    costs = None
    if cost:
        data = edge_data["cost"].values
        costs = sp.coo_matrix((data, (a0, a1)))

    cap = None
    if capacity:
        data = edge_data["capacity"].values
        cap = sp.coo_matrix((data, (a0, a1)))

    dem = None
    if demand:
        dem = node_data["demand"].values

    return G, cap, costs, dem
