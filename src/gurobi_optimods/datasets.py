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


def load_filepath(filename):
    file = str(DATA_FILE_DIR) + "/opf/" + filename
    return file


def load_caseopfmat(number):
    case = str(DATA_FILE_DIR) + "/opf/case" + number + ".mat"
    return case


def load_caseNYopf():  # real world data case
    case = str(DATA_FILE_DIR) + "/opf/caseNY.mat"
    return case


def load_opfdictcase():
    # Not quite equivalent to this?
    # from gurobi_optimods.opf import read_case_from_mat_file
    # return read_case_from_mat_file(load_caseopfmat("9"))
    casefile_dict = {
        "baseMVA": 100.0,
        "bus": {
            1: {
                "bus_i": 1,
                "type": 3,
                "Pd": 0.0,
                "Qd": 0.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "area": 0.0,
                "Vm": 1.0,
                "Va": 1.0,
                "baseKV": 345.0,
                "zone": 1.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
            },
            2: {
                "bus_i": 2,
                "type": 2,
                "Pd": 0.0,
                "Qd": 0.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "area": 0.0,
                "Vm": 1.0,
                "Va": 1.0,
                "baseKV": 345.0,
                "zone": 1.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
            },
            3: {
                "bus_i": 3,
                "type": 2,
                "Pd": 0.0,
                "Qd": 0.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "area": 0.0,
                "Vm": 1.0,
                "Va": 1.0,
                "baseKV": 345.0,
                "zone": 1.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
            },
            4: {
                "bus_i": 4,
                "type": 1,
                "Pd": 0.0,
                "Qd": 0.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "area": 0.0,
                "Vm": 1.0,
                "Va": 1.0,
                "baseKV": 345.0,
                "zone": 1.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
            },
            5: {
                "bus_i": 5,
                "type": 1,
                "Pd": 90.0,
                "Qd": 30.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "area": 0.0,
                "Vm": 1.0,
                "Va": 1.0,
                "baseKV": 345.0,
                "zone": 1.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
            },
            6: {
                "bus_i": 6,
                "type": 1,
                "Pd": 0.0,
                "Qd": 0.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "area": 0.0,
                "Vm": 1.0,
                "Va": 1.0,
                "baseKV": 345.0,
                "zone": 1.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
            },
            7: {
                "bus_i": 7,
                "type": 1,
                "Pd": 100.0,
                "Qd": 35.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "area": 0.0,
                "Vm": 1.0,
                "Va": 1.0,
                "baseKV": 345.0,
                "zone": 1.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
            },
            8: {
                "bus_i": 8,
                "type": 1,
                "Pd": 0.0,
                "Qd": 0.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "area": 0.0,
                "Vm": 1.0,
                "Va": 1.0,
                "baseKV": 345.0,
                "zone": 1.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
            },
            9: {
                "bus_i": 9,
                "type": 1,
                "Pd": 125.0,
                "Qd": 50.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "area": 0.0,
                "Vm": 1.0,
                "Va": 1.0,
                "baseKV": 345.0,
                "zone": 1.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
            },
        },
        "gen": {
            1: {
                "bus": 1,
                "Pg": 0.0,
                "Qg": 0.0,
                "Qmax": 300.0,
                "Qmin": -300.0,
                "Vg": 1,
                "mBase": 100,
                "status": 1,
                "Pmax": 250.0,
                "Pmin": 10.0,
                "Pc1": 0,
                "Pc2": 0,
                "Qc1min": 0,
                "Qc1max": 0,
                "Qc2min": 0,
                "Qc2max": 0,
                "ramp_agc": 0,
                "ramp_10": 0,
                "ramp_30": 0,
                "ramp_q": 0,
                "apf": 0,
            },
            2: {
                "bus": 2,
                "Pg": 163.0,
                "Qg": 0.0,
                "Qmax": 300.0,
                "Qmin": -300.0,
                "Vg": 1,
                "mBase": 100,
                "status": 1,
                "Pmax": 300.0,
                "Pmin": 10.0,
                "Pc1": 0,
                "Pc2": 0,
                "Qc1min": 0,
                "Qc1max": 0,
                "Qc2min": 0,
                "Qc2max": 0,
                "ramp_agc": 0,
                "ramp_10": 0,
                "ramp_30": 0,
                "ramp_q": 0,
                "apf": 0,
            },
            3: {
                "bus": 3,
                "Pg": 85.0,
                "Qg": 0.0,
                "Qmax": 300.0,
                "Qmin": -300.0,
                "Vg": 1,
                "mBase": 100,
                "status": 1,
                "Pmax": 270.0,
                "Pmin": 10.0,
                "Pc1": 0,
                "Pc2": 0,
                "Qc1min": 0,
                "Qc1max": 0,
                "Qc2min": 0,
                "Qc2max": 0,
                "ramp_agc": 0,
                "ramp_10": 0,
                "ramp_30": 0,
                "ramp_q": 0,
                "apf": 0,
            },
        },
        "branch": {
            1: {
                "fbus": 1,
                "tbus": 4,
                "r": 0.0,
                "x": 0.0576,
                "b": 0.0,
                "rateA": 250.0,
                "rateB": 250.0,
                "rateC": 250.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "angmin": -360.0,
                "angmax": 360.0,
            },
            2: {
                "fbus": 4,
                "tbus": 5,
                "r": 0.017,
                "x": 0.092,
                "b": 0.158,
                "rateA": 250.0,
                "rateB": 250.0,
                "rateC": 250.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "angmin": -360.0,
                "angmax": 360.0,
            },
            3: {
                "fbus": 5,
                "tbus": 6,
                "r": 0.039,
                "x": 0.17,
                "b": 0.358,
                "rateA": 150.0,
                "rateB": 150.0,
                "rateC": 150.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "angmin": -360.0,
                "angmax": 360.0,
            },
            4: {
                "fbus": 3,
                "tbus": 6,
                "r": 0.0,
                "x": 0.0586,
                "b": 0.0,
                "rateA": 300.0,
                "rateB": 300.0,
                "rateC": 300.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "angmin": -360.0,
                "angmax": 360.0,
            },
            5: {
                "fbus": 6,
                "tbus": 7,
                "r": 0.0119,
                "x": 0.1008,
                "b": 0.209,
                "rateA": 150.0,
                "rateB": 150.0,
                "rateC": 150.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "angmin": -360.0,
                "angmax": 360.0,
            },
            6: {
                "fbus": 7,
                "tbus": 8,
                "r": 0.0085,
                "x": 0.072,
                "b": 0.149,
                "rateA": 250.0,
                "rateB": 250.0,
                "rateC": 250.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "angmin": -360.0,
                "angmax": 360.0,
            },
            7: {
                "fbus": 8,
                "tbus": 2,
                "r": 0.0,
                "x": 0.0625,
                "b": 0.0,
                "rateA": 250.0,
                "rateB": 250.0,
                "rateC": 250.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "angmin": -360.0,
                "angmax": 360.0,
            },
            8: {
                "fbus": 8,
                "tbus": 9,
                "r": 0.032,
                "x": 0.161,
                "b": 0.306,
                "rateA": 250.0,
                "rateB": 250.0,
                "rateC": 250.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "angmin": -360.0,
                "angmax": 360.0,
            },
            9: {
                "fbus": 9,
                "tbus": 4,
                "r": 0.01,
                "x": 0.085,
                "b": 0.176,
                "rateA": 250.0,
                "rateB": 250.0,
                "rateC": 250.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "angmin": -360.0,
                "angmax": 360.0,
            },
        },
        "gencost": {
            1: {
                "costtype": 2,
                "startup": 1500,
                "shutdown": 0,
                "n": 3,
                "costvector": [0.11, 5, 150],
            },
            2: {
                "costtype": 2,
                "startup": 2000,
                "shutdown": 0,
                "n": 3,
                "costvector": [0.085, 1.2, 600],
            },
            3: {
                "costtype": 2,
                "startup": 3000,
                "shutdown": 0,
                "n": 3,
                "costvector": [0.1225, 1, 335],
            },
        },
    }
    return casefile_dict


def load_case9branchswitching():
    # we alter the original case dictionary in order to
    # create an artifical case, where turning off 2 branches
    # produces a better solution
    casefile_dict = load_opfdictcase()
    casefile_dict["branch"][10] = {
        "fbus": 1,
        "tbus": 2,
        "r": 0.0,
        "x": 0.0576,
        "b": 0.0,
        "rateA": 250.0,
        "rateB": 250.0,
        "rateC": 250.0,
        "ratio": 1.0,
        "angle": 0.0,
        "status": 1,
        "angmin": -360.0,
        "angmax": 360.0,
    }
    casefile_dict["branch"][11] = {
        "fbus": 1,
        "tbus": 3,
        "r": 0.0,
        "x": 0.0576,
        "b": 0.0,
        "rateA": 250.0,
        "rateB": 250.0,
        "rateC": 250.0,
        "ratio": 1.0,
        "angle": 0.0,
        "status": 1,
        "angmin": -360.0,
        "angmax": 360.0,
    }
    casefile_dict["branch"][12] = {
        "fbus": 2,
        "tbus": 3,
        "r": 0.0,
        "x": 0.0576,
        "b": 0.0,
        "rateA": 250.0,
        "rateB": 250.0,
        "rateC": 250.0,
        "ratio": 1.0,
        "angle": 0.0,
        "status": 1,
        "angmin": -360.0,
        "angmax": 360.0,
    }
    casefile_dict["gencost"][2]["costvector"] = [0.85, 10.2, 1200]

    return casefile_dict
