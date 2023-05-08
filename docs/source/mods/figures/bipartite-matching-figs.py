from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from gurobi_optimods.bipartite_matching import maximum_bipartite_matching

here = Path(__file__).parent


def matching_example():
    # side by side graph with its maximum matching

    n = 3
    m = 4

    graph = nx.Graph()
    graph.add_nodes_from(list(range(n + m)))
    graph.add_edges_from([(0, 3), (0, 4), (0, 6), (1, 3), (2, 4), (2, 5)])
    nodes1 = np.arange(3)
    nodes2 = np.arange(3, 3 + 4)

    matching = maximum_bipartite_matching(graph, nodes1, nodes2)
    fig, (ax1, ax2) = plt.subplots(1, 2)
    layout = {
        0: [-1, -0.75],
        1: [-1, 0],
        2: [-1, 0.75],
        3: [1, -0.75],
        4: [1, -0.25],
        5: [1, 0.25],
        6: [1, 0.75],
    }
    nx.draw(graph, layout, ax=ax1)
    nx.draw(matching, layout, ax=ax2)
    fig.savefig(here / "bipartite-matching-example.png")


def matching_network():
    # equivalent max flow network to solve

    n = 3
    m = 4
    graph = nx.DiGraph()
    graph.add_nodes_from(list(range(n + m)))
    graph.add_edges_from([(0, 3), (0, 4), (0, 6), (1, 3), (2, 4), (2, 5)])
    source = n + m
    sink = n + m + 1
    graph.add_nodes_from([source, sink])
    for i in range(n):
        graph.add_edge(source, i)
    for j in range(n, n + m):
        graph.add_edge(j, sink)

    layout = {
        0: [-1, -0.75],
        1: [-1, 0],
        2: [-1, 0.75],
        3: [1, -0.75],
        4: [1, -0.25],
        5: [1, 0.25],
        6: [1, 0.75],
        source: [-3.0, 0.0],
        sink: [3.0, 0.0],
    }
    fig = plt.figure()
    nx.draw(graph, layout, ax=plt.gca())
    fig.savefig(here / "bipartite-matching-flow.png")


if __name__ == "__main__":
    matching_example()
    matching_network()
