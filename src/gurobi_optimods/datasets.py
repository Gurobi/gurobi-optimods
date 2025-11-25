"""
Module for loading datasets for use in optimods examples, in the same vein
as sklearn.datasets.
"""

import json
import pathlib

import numpy as np
import pandas as pd
import scipy.sparse as sp

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


def load_siouxfalls_network_data():
    # edge and node data to create a graph
    edge_data = pd.read_csv(DATA_FILE_DIR / "graphs/siouxfalls_edges.csv")
    node_data = pd.read_csv(DATA_FILE_DIR / "graphs/siouxfalls_nodes.csv")

    # line data including, line-path, capacity, fixed cost, and operating cost
    linepath_data = pd.read_csv(DATA_FILE_DIR / "graphs/siouxfalls_linepaths.csv")
    line_data = pd.read_csv(DATA_FILE_DIR / "graphs/siouxfalls_lines.csv")

    # demand data
    demand_data = pd.read_csv(DATA_FILE_DIR / "graphs/siouxfalls_demand.csv")

    return (node_data, edge_data, line_data, linepath_data, demand_data)


def load_berlin_metro_reduced_graph_data():
    from networkx.readwrite import json_graph

    # read graph
    with open(DATA_FILE_DIR / "graphs/uberlin_reduced_graph.json", "r") as f:
        data = json.load(f)

    # Convert the JSON data to a NetworkX graph
    graph = json_graph.node_link_graph(data, edges="links")

    # line path data
    linepath_data = pd.read_csv(DATA_FILE_DIR / "graphs/uberlin_reduced_linepaths.csv")

    return (graph, linepath_data)


def load_berlin_metro_graph_data():
    from networkx.readwrite import json_graph

    # read graph
    with open(DATA_FILE_DIR / "graphs/uberlin_graph.json", "r") as f:
        data = json.load(f)

    # Convert the JSON data to a NetworkX graph
    graph = json_graph.node_link_graph(data, edges="links")

    # line path data
    linepath_data = pd.read_csv(DATA_FILE_DIR / "graphs/uberlin_linepaths.csv")

    return (graph, linepath_data)


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
    edge_data,
    node_data,
    capacity=True,
    cost=True,
    demand=True,
    use_multigraph=False,
):
    """
    Convert from a pandas DataFrame to a networkx.MultiDiGraph with the appropriate
    attributes. For edges: `capacity`, and `cost`. For nodes: `demand`.
    """
    import networkx as nx

    graph_type = nx.MultiDiGraph if use_multigraph else nx.DiGraph

    G = nx.from_pandas_edgelist(
        edge_data.reset_index(), create_using=graph_type(), edge_attr=True
    )
    if demand:
        for i, d in node_data.iterrows():
            G.add_node(i, demand=d.demand)
    return G


def _convert_pandas_to_scipy(
    edge_data, node_data, capacity=True, cost=True, demand=True
):
    """
    Convert from a pandas DataFrame to several scipy.sparse.coo_array contain
    the graph structure, the capacity and cost values per edge, and the demand
    values per node.
    """
    coords = edge_data.index.to_numpy()

    a0 = np.array([c[0] for c in coords])
    a1 = np.array([c[1] for c in coords])

    data = np.ones(len(coords), dtype=np.int64)
    G = sp.coo_array((data, (a0, a1)))

    costs = None
    if cost:
        data = edge_data["cost"].values
        costs = sp.coo_array((data, (a0, a1)))

    cap = None
    if capacity:
        data = edge_data["capacity"].values
        cap = sp.coo_array((data, (a0, a1)))

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


def load_facility_location():
    """
    Load facility location example dataset.

    This dataset contains a realistic example with 15 customers and 8 potential
    facility locations across 4 US regions. The data is based on real US city
    coordinates with transportation costs computed using great-circle distances.

    Returns
    -------
    AttrDict
        Dictionary-like object with the following keys:

        - ``customer_data``: DataFrame with columns ``customer`` and ``demand``
        - ``facility_data``: DataFrame with columns ``facility``, ``capacity``,
          ``fixed_cost``, and ``region`` (for fairness constraints)
        - ``transportation_cost``: DataFrame with columns ``customer``, ``facility``,
          and ``cost`` (computed from geographic distances)

    Examples
    --------
    >>> from gurobi_optimods import datasets
    >>> data = datasets.load_facility_location()
    >>> data.customer_data.shape
    (15, 2)
    >>> data.facility_data.shape
    (8, 4)
    >>> 'region' in data.facility_data.columns
    True
    """
    return AttrDict(
        customer_data=pd.read_csv(
            DATA_FILE_DIR / "facility_location/customer_data.csv"
        ),
        facility_data=pd.read_csv(
            DATA_FILE_DIR / "facility_location/facility_data.csv"
        ),
        transportation_cost=pd.read_csv(
            DATA_FILE_DIR / "facility_location/transportation_cost.csv"
        ),
    )
