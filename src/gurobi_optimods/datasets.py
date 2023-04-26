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
            DATA_FILE_DIR / "workforce/availability.csv", parse_dates=["Shift"]
        ),
        pay_rates=pd.read_csv(DATA_FILE_DIR / "workforce/pay_rates.csv"),
        shift_requirements=pd.read_csv(
            DATA_FILE_DIR / "workforce/shift_requirements.csv", parse_dates=["Shift"]
        ),
    )


def _convert_pandas_to_digraph(edge_data, node_data):
    """
    Convert from a pandas DataFrame to a networkx.DiGraph with the appropriate
    attributes. For edges: `capacity`, and `cost`. For nodes: `demand`.
    """
    G = nx.from_pandas_edgelist(
        edge_data.reset_index(), create_using=nx.DiGraph(), edge_attr=True
    )
    for i, d in node_data.iterrows():
        G.add_node(i, demand=d.demand)
    return G


def _convert_pandas_to_scipy(edge_data, node_data):
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

    data = edge_data["capacity"].values
    cap = sp.coo_matrix((data, (a0, a1)))

    data = edge_data["cost"].values
    costs = sp.coo_matrix((data, (a0, a1)))

    return G, cap, costs, node_data["demand"].values


def load_min_cost_flow(drop_pos=True):
    edge_data = pd.read_csv(DATA_FILE_DIR / "graphs/edge_data1.csv").set_index(
        ["source", "target"]
    )
    node_data = pd.read_csv(DATA_FILE_DIR / "graphs/node_data1.csv", index_col=0)
    if drop_pos:
        node_data.drop(columns=["posx", "posy"], inplace=True)
    return edge_data, node_data


def load_min_cost_flow_networkx():
    edge_data, node_data = load_min_cost_flow()
    return _convert_pandas_to_digraph(edge_data, node_data)


def load_min_cost_flow_scipy():
    edge_data, node_data = load_min_cost_flow()
    return _convert_pandas_to_scipy(edge_data, node_data)


def load_min_cost_flow2():
    return (
        pd.read_csv(DATA_FILE_DIR / "graphs/edge_data2.csv").set_index(
            ["source", "target"]
        ),
        pd.read_csv(DATA_FILE_DIR / "graphs/node_data2.csv", index_col=0),
    )


def load_min_cost_flow2_networkx():
    edge_data, node_data = load_min_cost_flow2()
    return _convert_pandas_to_digraph(edge_data, node_data)


def load_min_cost_flow2_scipy():
    edge_data, node_data = load_min_cost_flow2()
    return _convert_pandas_to_scipy(edge_data, node_data)


def load_diet():
    return AttrDict(
        categories=pd.read_csv(DATA_FILE_DIR / "diet-categories.csv"),
        foods=pd.read_csv(DATA_FILE_DIR / "diet-foods.csv"),
        nutrition_values=pd.read_csv(DATA_FILE_DIR / "diet-values.csv"),
    )
