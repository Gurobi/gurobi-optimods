"""
Module for loading datasets for use in optimods examples, in the same vein
as sklearn.datasets.
"""

import csv
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


def _load_simple_graph_pandas(drop_pos=True, capacity=True, cost=True, demand=True):
    edge_data = pd.read_csv(DATA_FILE_DIR / "graphs/simple_graph_edges.csv").set_index(
        ["source", "target"]
    )
    node_data = pd.read_csv(
        DATA_FILE_DIR / "graphs/simple_graph_nodes.csv", index_col=0
    )
    if drop_pos:
        node_data.drop(columns=["posx", "posy"], inplace=True)
    if not capacity:
        edge_data.drop(columns=["capacity"], inplace=True)
    if not cost:
        edge_data.drop(columns=["cost"], inplace=True)
    if not demand:
        node_data.drop(columns=["demand"], inplace=True)
    return edge_data, node_data


def simple_graph_pandas():
    return _load_simple_graph_pandas(capacity=True, cost=True, demand=True)


def simple_graph_networkx():
    edge_data, node_data = _load_simple_graph_pandas(
        capacity=True, cost=True, demand=True
    )
    return _convert_pandas_to_digraph(
        edge_data, node_data, capacity=True, cost=True, demand=True
    )


def simple_graph_scipy(capacity=True, cost=True, demand=True):
    edge_data, node_data = _load_simple_graph_pandas(
        capacity=True, cost=True, demand=True
    )
    return _convert_pandas_to_scipy(
        edge_data, node_data, capacity=True, cost=True, demand=True
    )


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


def load_sharpe_ratio():
    data = pd.read_csv(DATA_FILE_DIR / "sharpe-ratio/log-returns.csv", index_col=0)
    # Annualize covariance-variance matrix and expected returns
    return AttrDict(cov_matrix=data.cov() * len(data.index), mu=data.sum())


def load_opf_example(case):
    from gurobi_optimods.opf.io import read_case_matpower

    file_path = DATA_FILE_DIR / f"opf/{case}.mat"
    return read_case_matpower(file_path)


def load_opf_extra(extra):
    data_files = {
        "case9-coordinates": ("case9coords.csv", "coords"),
        "case9-voltages": ("case9volts.csv", "volts"),
        "caseNY-coordinates": ("nybuses.csv", "coords"),
    }

    file_name, file_type = data_files[extra]
    file_path = DATA_FILE_DIR.joinpath("opf").joinpath(file_name)
    data = pd.read_csv(file_path)

    if file_type == "coords":

        def mapper(row):
            return row["bus_i"], (row["lat"], row["lon"])

    elif file_type == "volts":

        def mapper(row):
            return row["bus_i"], (row["Vm"], row["Va"])

    return dict(mapper(record) for record in data.to_dict("records"))
