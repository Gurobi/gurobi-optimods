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
    )


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
