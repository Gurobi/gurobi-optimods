"""
Module for loading datasets for use in optimods examples, in the same vein
as sklearn.datasets.
"""

import pathlib

import pandas as pd

try:
    import networkx as nx
except ImportError:
    nx = None

from gurobi_optimods.network_util import (
    _convert_pandas_to_digraph,
    _convert_pandas_to_scipy,
)

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
        preferences=pd.read_csv(
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


def load_diet():
    return AttrDict(
        categories=pd.read_csv(DATA_FILE_DIR / "diet-categories.csv"),
        foods=pd.read_csv(DATA_FILE_DIR / "diet-foods.csv"),
        nutrition_values=pd.read_csv(DATA_FILE_DIR / "diet-values.csv"),
    )
