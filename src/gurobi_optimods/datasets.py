"""
Module for loading datasets for use in optimods examples, in the same vein
as sklearn.datasets.
"""

import pathlib

import pandas as pd
import networkx as nx
import matplotlib
import matplotlib.pyplot

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
            DATA_FILE_DIR / "availability.csv", parse_dates=["Shift"]
        ),
        pay_rates=pd.read_csv(DATA_FILE_DIR / "pay_rates.csv"),
        shift_requirements=pd.read_csv(
            DATA_FILE_DIR / "shift_requirements.csv", parse_dates=["Shift"]
        ),
    )


def load_network_flow_example_data():
    return nx.read_gml(DATA_FILE_DIR / "network_flow.gml", destringizer=int), 0, 4


def load_network_flow_example_data2():
    return nx.read_gml(DATA_FILE_DIR / "network_flow2.gml", destringizer=int), "s", "t"


def _create_network_flow2():
    G = nx.DiGraph()

    G.add_node("s", pos=(0, 0), demand=2)
    G.add_node(1, pos=(1, 0.5), demand=0)
    G.add_node(2, pos=(1, -0.5), demand=1)
    G.add_node(3, pos=(2, 0.5), demand=-1)
    G.add_node(4, pos=(2, -0.5), demand=0)
    G.add_node("t", pos=(3, 0), demand=-2)

    G.add_edge("s", 1, capacity=2, cost=9)
    G.add_edge("s", 2, capacity=2, cost=7)
    G.add_edge(1, 3, capacity=1, cost=1)
    G.add_edge(2, 3, capacity=1, cost=10)
    G.add_edge(2, 4, capacity=2, cost=6)
    G.add_edge(3, "t", capacity=2, cost=1)
    G.add_edge(4, "t", capacity=2, cost=1)

    # TODO: plot and place in docs
    fig = matplotlib.pyplot.figure()
    pos = nx.get_node_attributes(G, "pos")
    lab = {}
    for e in G.edges():
        lab[e] = str((G.edges[e]["capacity"], G.edges[e]["cost"]))
    print(lab)
    nx.draw(G, ax=fig.add_subplot(), with_labels=True, pos=pos)
    nx.draw_networkx_edge_labels(G, pos, lab)
    # Save plot to file
    matplotlib.use("Agg")
    fig.savefig("graph.png", dpi=1000)

    nx.write_gml(G, "network_flow.gml")


def _create_network_flow1():
    G = nx.DiGraph()

    G.add_edge(0, 1, capacity=15, cost=4)
    G.add_edge(0, 2, capacity=8, cost=4)
    G.add_edge(1, 3, capacity=4, cost=2)
    G.add_edge(1, 2, capacity=20, cost=2)
    G.add_edge(1, 4, capacity=10, cost=6)
    G.add_edge(2, 3, capacity=15, cost=1)
    G.add_edge(3, 4, capacity=20, cost=2)
    G.add_edge(2, 4, capacity=5, cost=3)
    G.add_edge(4, 2, capacity=4, cost=3)

    nx.set_node_attributes(G, 0, "demand")

    G.nodes[0]["demand"] = 20
    G.nodes[3]["demand"] = -5
    G.nodes[4]["demand"] = -15

    # TODO: plot and place in docs
    fig = matplotlib.pyplot.figure()
    nx.draw(G, ax=fig.add_subplot())
    # Save plot to file
    matplotlib.use("Agg")
    fig.savefig("graph.png")

    nx.write_gml(G, "network_flow.gml")
