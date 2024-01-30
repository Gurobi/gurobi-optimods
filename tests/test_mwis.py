import unittest
from itertools import combinations

import numpy as np
import pandas as pd
import scipy.sparse as sp
from numpy.testing import assert_array_equal

try:
    import networkx as nx
except ImportError:
    nx = None

from gurobi_optimods.mwis import (
    maximum_weighted_clique,
    maximum_weighted_independent_set,
)


def get_graph(num_vertices, density, seed, graph_type, scipy_func=sp.csr_matrix):
    """Create the upper triangular adjacency matrix for a random graph
        with the given number of vertices and density.

    Args:
        num_vertices (int): Number of vertices.
        density (float): Graph density between 0 and 1.
        seed (int): Random seed value.
        graph_type (spmatrix, Graph, or DataFrame): Graph type to return.
        scipy_func (func): Sparse func to be used if graph_type is spmatrix.
            Defaults to sp.csr_matrix

    Returns:
        graph: The spmatrix, Graph, or DataFrame graph.
    """
    rng = np.random.default_rng(seed=seed)
    num_edges = int(0.5 * density * num_vertices * (num_vertices - 1))
    if num_edges == 0:
        if graph_type == "spmatrix":
            return scipy_func((num_vertices, num_vertices))
        elif graph_type == "Graph":
            graph = nx.Graph()
            graph.add_nodes_from(range(num_vertices))
            return graph
        elif graph_type == "DataFrame":
            return pd.DataFrame({}, columns=["node1", "node2"])

    edges = rng.choice(
        list(combinations(range(num_vertices), 2)), num_edges, replace=False
    )
    if graph_type == "spmatrix":
        rows, cols = zip(*edges)
        data = np.ones(num_edges)
        return scipy_func((data, (rows, cols)), shape=(num_vertices, num_vertices))
    elif graph_type == "Graph":
        graph = nx.Graph()
        graph.add_nodes_from(range(num_vertices))
        graph.add_edges_from(edges)
        return graph
    elif graph_type == "DataFrame":
        return pd.DataFrame(edges, columns=["node1", "node2"])


class TestMWISScipySparse(unittest.TestCase):
    def test_random_graph(self):
        for density in [np.random.random() for _ in range(5)]:
            for func in [
                sp.csr_matrix,
                sp.csc_matrix,
                sp.coo_matrix,
                sp.csr_array,
                sp.csc_array,
                sp.coo_array,
            ]:
                with self.subTest(density=density, func=func):
                    num_vertices, density, seed = 10, density, 0
                    graph = get_graph(num_vertices, density, seed, "spmatrix", func)
                    weights = np.random.randint(1, 100, size=num_vertices)
                    mwis = maximum_weighted_independent_set(graph, weights)
                    self.assertLessEqual(len(mwis.x), num_vertices)
                    self.assertGreaterEqual(len(mwis.x), 1)
                    self.assertLessEqual(mwis.f, weights.sum())
                    self.assertGreaterEqual(mwis.f, weights.max())

    def test_complete_graph(self):
        for num_vertices in [np.random.randint(10, 20) for _ in range(5)]:
            density, seed = 1, 0
            with self.subTest(num_vertices=num_vertices):
                graph = get_graph(num_vertices, density, seed, "spmatrix")
                weights = np.random.randint(1, 100, size=num_vertices)
                mwis = maximum_weighted_independent_set(graph, weights)
                self.assertEqual(len(mwis.x), 1)
                self.assertEqual(mwis.f, weights.max())

    def test_empty_graph(self):
        for num_vertices in [np.random.randint(10, 20) for _ in range(5)]:
            density, seed = 0, 0
            with self.subTest(num_vertices=num_vertices):
                graph = get_graph(num_vertices, density, seed, "spmatrix")
                weights = np.random.randint(1, 100, size=num_vertices)
                mwis = maximum_weighted_independent_set(graph, weights)
                self.assertEqual(len(mwis.x), num_vertices)
                self.assertEqual(mwis.f, weights.sum())

    def test_known_graph(self):
        rows = [0, 0, 0, 1, 1, 2, 2, 3, 4, 4, 5, 6]
        cols = [1, 3, 4, 3, 5, 3, 6, 7, 5, 7, 6, 7]
        data = [1 for _ in range(12)]
        graph = sp.csr_matrix((data, (rows, cols)), shape=(8, 8))
        weights = np.array([2**i for i in range(8)])
        mwis = maximum_weighted_independent_set(graph, weights)
        assert_array_equal(mwis.x, np.array([0, 2, 5, 7]))
        self.assertEqual(mwis.f, 165)

    def test_solver_params(self):
        num_vertices, density, seed = 10, 0.75, 0
        graph = get_graph(num_vertices, density, seed, "spmatrix")
        weights = np.random.randint(1, 100, size=num_vertices)
        mwis = maximum_weighted_independent_set(
            graph, weights, solver_params={"presolve": 2, "method": 1}
        )
        self.assertGreaterEqual(len(mwis.x), 1)
        self.assertLessEqual(mwis.f, weights.sum())


class TestMWCScipySparse(unittest.TestCase):
    def test_random_graph(self):
        for density in [np.random.random() for _ in range(5)]:
            for func in [
                sp.csr_matrix,
                sp.csc_matrix,
                sp.coo_matrix,
                sp.csr_array,
                sp.csc_array,
                sp.coo_array,
            ]:
                with self.subTest(density=density, func=func):
                    num_vertices, density, seed = 10, density, 0
                    graph = get_graph(num_vertices, density, seed, "spmatrix", func)
                    weights = np.random.randint(1, 100, size=num_vertices)
                    mwc = maximum_weighted_clique(graph, weights)
                    self.assertLessEqual(len(mwc.x), num_vertices)
                    self.assertGreaterEqual(len(mwc.x), 1)
                    self.assertLessEqual(mwc.f, weights.sum())
                    self.assertGreaterEqual(mwc.f, weights.max())

    def test_complete_graph(self):
        for num_vertices in [np.random.randint(10, 20) for _ in range(5)]:
            density, seed = 1, 0
            with self.subTest(num_vertices=num_vertices):
                graph = get_graph(num_vertices, density, seed, "spmatrix")
                weights = np.random.randint(1, 100, size=num_vertices)
                mwc = maximum_weighted_clique(graph, weights)
                self.assertEqual(len(mwc.x), num_vertices)
                self.assertEqual(mwc.f, weights.sum())

    def test_empty_graph(self):
        for num_vertices in [np.random.randint(10, 20) for _ in range(5)]:
            density, seed = 0, 0
            with self.subTest(num_vertices=num_vertices):
                graph = get_graph(num_vertices, density, seed, "spmatrix")
                weights = np.random.randint(1, 100, size=num_vertices)
                mwc = maximum_weighted_clique(graph, weights)
                self.assertEqual(len(mwc.x), 1)
                self.assertEqual(mwc.f, weights.max())

    def test_known_graph(self):
        rows = [0, 0, 0, 1, 1, 2, 2, 3, 4, 4, 5, 6]
        cols = [1, 3, 4, 3, 5, 3, 6, 7, 5, 7, 6, 7]
        data = [1 for _ in range(12)]
        adjacency_matrix = sp.csr_matrix((data, (rows, cols)), shape=(8, 8))
        weights = np.array([2**i for i in range(8)])
        mwc = maximum_weighted_clique(adjacency_matrix, weights)
        assert_array_equal(mwc.x, np.array([6, 7]))
        self.assertEqual(mwc.f, 192)

    def test_solver_params(self):
        num_vertices, density, seed = 10, 0.75, 0
        graph = get_graph(num_vertices, density, seed, "spmatrix")
        weights = np.random.randint(1, 100, size=num_vertices)
        mwc = maximum_weighted_clique(
            graph, weights, solver_params={"presolve": 2, "method": 1}
        )
        self.assertGreaterEqual(len(mwc.x), 1)
        self.assertLessEqual(mwc.f, weights.sum())


class TestMISPandas(unittest.TestCase):
    def test_random_graph(self):
        for density in [np.random.random() for _ in range(5)]:
            for seed in range(5):
                with self.subTest(density=density, seed=seed):
                    num_vertices, density, seed = 10, density, seed
                    graph = get_graph(num_vertices, density, seed, "DataFrame")
                    weights = pd.DataFrame(
                        np.random.randint(1, 100, size=num_vertices),
                        columns=["weights"],
                    )
                    mwis = maximum_weighted_independent_set(graph, weights)
                    self.assertLessEqual(len(mwis.x), num_vertices)
                    self.assertGreaterEqual(len(mwis.x), 1)
                    self.assertLessEqual(mwis.f, weights["weights"].sum())
                    self.assertGreaterEqual(mwis.f, weights["weights"].max())

    def test_complete_graph(self):
        for num_vertices in [np.random.randint(10, 20) for _ in range(5)]:
            for seed in range(5):
                density = 1
                with self.subTest(num_vertices=num_vertices, seed=seed):
                    graph = get_graph(num_vertices, density, seed, "DataFrame")
                    weights = pd.DataFrame(
                        np.random.randint(1, 100, size=num_vertices),
                        columns=["weights"],
                    )
                    mwis = maximum_weighted_independent_set(graph, weights)
                    self.assertEqual(len(mwis.x), 1)
                    self.assertEqual(mwis.f, weights["weights"].max())

    def test_empty_graph(self):
        for num_vertices in [np.random.randint(10, 20) for _ in range(5)]:
            for seed in range(5):
                density = 0
                with self.subTest(num_vertices=num_vertices, seed=seed):
                    graph = get_graph(num_vertices, density, seed, "DataFrame")
                    weights = pd.DataFrame(
                        np.random.randint(1, 100, size=num_vertices),
                        columns=["weights"],
                    )
                    mwis = maximum_weighted_independent_set(graph, weights)
                    self.assertEqual(len(mwis.x), num_vertices)
                    self.assertEqual(mwis.f, weights["weights"].sum())

    def test_known_graph(self):
        edges = [
            (0, 1),
            (0, 3),
            (0, 4),
            (1, 3),
            (1, 5),
            (2, 3),
            (2, 6),
            (3, 7),
            (4, 5),
            (4, 7),
            (5, 6),
            (6, 7),
        ]
        frame = pd.DataFrame(edges, columns=["node1", "node2"])
        weights = pd.DataFrame(
            np.array([2**i for i in range(8)]), columns=["weights"]
        )
        mwis = maximum_weighted_independent_set(frame, weights)
        assert_array_equal(mwis.x, np.array([0, 2, 5, 7]))
        self.assertEqual(mwis.f, 165)

    def test_solver_params(self):
        num_vertices, density, seed = 10, 0.75, 0
        frame = get_graph(num_vertices, density, seed, "DataFrame")
        weights = pd.DataFrame(
            np.random.randint(1, 100, size=num_vertices),
            columns=["weights"],
        )
        mwis = maximum_weighted_independent_set(
            frame, weights, solver_params={"presolve": 2, "method": 1}
        )
        self.assertGreaterEqual(len(mwis.x), 1)
        self.assertLessEqual(mwis.f, weights["weights"].sum())


class TestMWCPandas(unittest.TestCase):
    def test_random_graph(self):
        for density in [np.random.random() for _ in range(5)]:
            for seed in range(5):
                with self.subTest(density=density, seed=seed):
                    num_vertices, density, seed = 10, density, seed
                    graph = get_graph(num_vertices, density, seed, "DataFrame")
                    weights = pd.DataFrame(
                        np.random.randint(1, 100, size=num_vertices),
                        columns=["weights"],
                    )
                    mwc = maximum_weighted_clique(graph, weights)
                    self.assertLessEqual(len(mwc.x), num_vertices)
                    self.assertGreaterEqual(len(mwc.x), 1)
                    self.assertLessEqual(mwc.f, weights["weights"].sum())
                    self.assertGreaterEqual(mwc.f, weights["weights"].max())

    def test_complete_graph(self):
        for num_vertices in [np.random.randint(10, 20) for _ in range(5)]:
            for seed in range(5):
                density = 1
                with self.subTest(num_vertices=num_vertices, seed=seed):
                    graph = get_graph(num_vertices, density, seed, "DataFrame")
                    weights = pd.DataFrame(
                        np.random.randint(1, 100, size=num_vertices),
                        columns=["weights"],
                    )
                    mwc = maximum_weighted_clique(graph, weights)
                    self.assertEqual(len(mwc.x), num_vertices)
                    self.assertEqual(mwc.f, weights["weights"].sum())

    def test_empty_graph(self):
        for num_vertices in [np.random.randint(10, 20) for _ in range(5)]:
            for seed in range(5):
                density = 0
                with self.subTest(num_vertices=num_vertices, seed=seed):
                    graph = get_graph(num_vertices, density, seed, "DataFrame")
                    weights = pd.DataFrame(
                        np.random.randint(1, 100, size=num_vertices),
                        columns=["weights"],
                    )
                    mwc = maximum_weighted_clique(graph, weights)
                    self.assertEqual(len(mwc.x), 1)
                    self.assertEqual(mwc.f, weights["weights"].max())

    def test_known_graph(self):
        edges = [
            (0, 1),
            (0, 3),
            (0, 4),
            (1, 3),
            (1, 5),
            (2, 3),
            (2, 6),
            (3, 7),
            (4, 5),
            (4, 7),
            (5, 6),
            (6, 7),
        ]
        frame = pd.DataFrame(edges, columns=["node1", "node2"])
        weights = pd.DataFrame(
            np.array([2**i for i in range(8)]), columns=["weights"]
        )
        mwc = maximum_weighted_clique(frame, weights)
        assert_array_equal(mwc.x, np.array([6, 7]))
        self.assertEqual(mwc.f, 192)

    def test_solver_params(self):
        num_vertices, density, seed = 10, 0.75, 0
        frame = get_graph(num_vertices, density, seed, "DataFrame")
        weights = pd.DataFrame(
            np.random.randint(1, 100, size=num_vertices),
            columns=["weights"],
        )
        mwc = maximum_weighted_independent_set(
            frame, weights, solver_params={"presolve": 2, "method": 1}
        )
        self.assertGreaterEqual(len(mwc.x), 1)
        self.assertLessEqual(mwc.f, weights["weights"].sum())


@unittest.skipIf(nx is None, "networkx is not installed")
class TestMWISNetworkx(unittest.TestCase):
    def test_random_graph(self):
        for density in [np.random.random() for _ in range(5)]:
            for seed in range(5):
                with self.subTest(density=density, seed=seed):
                    num_vertices, density, seed = 10, density, seed
                    graph = get_graph(num_vertices, density, seed, "Graph")
                    weights = np.random.randint(1, 100, size=num_vertices)
                    mwis = maximum_weighted_independent_set(graph, weights)
                    self.assertLessEqual(len(mwis.x), num_vertices)
                    self.assertGreaterEqual(len(mwis.x), 1)
                    self.assertLessEqual(mwis.f, weights.sum())
                    self.assertGreaterEqual(mwis.f, weights.max())

    def test_complete_graph(self):
        for num_vertices in [np.random.randint(10, 20) for _ in range(5)]:
            for seed in range(5):
                density = 1
                with self.subTest(num_vertices=num_vertices, seed=seed):
                    graph = get_graph(num_vertices, density, seed, "Graph")
                    weights = np.random.randint(1, 100, size=num_vertices)
                    mwis = maximum_weighted_independent_set(graph, weights)
                    self.assertEqual(len(mwis.x), 1)
                    self.assertEqual(mwis.f, weights.max())

    def test_empty_graph(self):
        for num_vertices in [np.random.randint(10, 20) for _ in range(5)]:
            for seed in range(5):
                density = 0
                with self.subTest(num_vertices=num_vertices, seed=seed):
                    graph = get_graph(num_vertices, density, seed, "Graph")
                    weights = np.random.randint(1, 100, size=num_vertices)
                    mwis = maximum_weighted_independent_set(graph, weights)
                    self.assertEqual(len(mwis.x), num_vertices)
                    self.assertEqual(mwis.f, weights.sum())

    def test_known_graph(self):
        edges = [
            (0, 1),
            (0, 3),
            (0, 4),
            (1, 3),
            (1, 5),
            (2, 3),
            (2, 6),
            (3, 7),
            (4, 5),
            (4, 7),
            (5, 6),
            (6, 7),
        ]
        graph = nx.Graph()
        graph.add_edges_from(edges)
        weights = np.array([2**i for i in range(8)])
        mwis = maximum_weighted_independent_set(graph, weights)
        assert_array_equal(mwis.x, np.array([0, 2, 5, 7]))
        self.assertEqual(mwis.f, 165)

    def test_solver_params(self):
        num_vertices, density, seed = 10, 0.75, 0
        graph = get_graph(num_vertices, density, seed, "Graph")
        weights = np.random.randint(1, 100, size=num_vertices)
        mwis = maximum_weighted_independent_set(
            graph, weights, solver_params={"presolve": 2, "method": 1}
        )
        self.assertGreaterEqual(len(mwis.x), 1)
        self.assertLessEqual(mwis.f, weights.sum())


@unittest.skipIf(nx is None, "networkx is not installed")
class TestMWCNetworkx(unittest.TestCase):
    def test_random_graph(self):
        for density in [np.random.random() for _ in range(5)]:
            for seed in range(5):
                with self.subTest(density=density, seed=seed):
                    num_vertices, density, seed = 10, density, seed
                    graph = get_graph(num_vertices, density, seed, "Graph")
                    weights = np.random.randint(1, 100, size=num_vertices)
                    mwc = maximum_weighted_clique(graph, weights)
                    self.assertLessEqual(len(mwc.x), num_vertices)
                    self.assertGreaterEqual(len(mwc.x), 1)
                    self.assertLessEqual(mwc.f, weights.sum())
                    self.assertGreaterEqual(mwc.f, weights.max())

    def test_complete_graph(self):
        for num_vertices in [np.random.randint(10, 20) for _ in range(5)]:
            for seed in range(5):
                density = 1
                with self.subTest(num_vertices=num_vertices, seed=seed):
                    graph = get_graph(num_vertices, density, seed, "Graph")
                    weights = np.random.randint(1, 100, size=num_vertices)
                    mwc = maximum_weighted_clique(graph, weights)
                    self.assertEqual(len(mwc.x), num_vertices)
                    self.assertEqual(mwc.f, weights.sum())

    def test_empty_graph(self):
        for num_vertices in [np.random.randint(10, 20) for _ in range(5)]:
            for seed in range(5):
                density = 0
                with self.subTest(num_vertices=num_vertices, seed=seed):
                    graph = get_graph(num_vertices, density, seed, "Graph")
                    weights = np.random.randint(1, 100, size=num_vertices)
                    mwc = maximum_weighted_clique(graph, weights)
                    self.assertEqual(len(mwc.x), 1)
                    self.assertEqual(mwc.f, weights.max())

    def test_known_graph(self):
        edges = [
            (0, 1),
            (0, 3),
            (0, 4),
            (1, 3),
            (1, 5),
            (2, 3),
            (2, 6),
            (3, 7),
            (4, 5),
            (4, 7),
            (5, 6),
            (6, 7),
        ]
        graph = nx.Graph()
        graph.add_edges_from(edges)
        weights = np.array([2**i for i in range(8)])
        mwc = maximum_weighted_clique(graph, weights)
        assert_array_equal(mwc.x, np.array([6, 7]))
        self.assertEqual(mwc.f, 192)

    def test_solver_params(self):
        num_vertices, density, seed = 10, 0.75, 0
        graph = get_graph(num_vertices, density, seed, "Graph")
        weights = np.random.randint(1, 100, size=num_vertices)
        mwc = maximum_weighted_clique(
            graph, weights, solver_params={"presolve": 2, "method": 1}
        )
        self.assertGreaterEqual(len(mwc.x), 1)
        self.assertLessEqual(mwc.f, weights.sum())


class TestMWISAll(unittest.TestCase):
    def test_three_approaches(self):
        for density in [np.random.random() for _ in range(5)]:
            for num_vertices in [np.random.randint(0, 30) for _ in range(5)]:
                for seed in range(5):
                    with self.subTest(
                        num_vertices=num_vertices, density=density, seed=seed
                    ):
                        num_vertices, density, seed = 10, density, seed
                        graph_sp = get_graph(num_vertices, density, seed, "spmatrix")
                        graph_pd = get_graph(num_vertices, density, seed, "DataFrame")
                        weights = np.random.randint(1, 100, size=num_vertices)
                        weights_pd = pd.DataFrame(weights, columns=["weights"])

                        mwis_sp = maximum_weighted_independent_set(graph_sp, weights)
                        mwis_pd = maximum_weighted_independent_set(graph_pd, weights_pd)
                        self.assertEqual(mwis_sp.f, mwis_pd.f)

                        if nx is not None:
                            graph_nx = get_graph(num_vertices, density, seed, "Graph")
                            mwis_nx = maximum_weighted_independent_set(
                                graph_nx, weights
                            )
                            self.assertEqual(mwis_sp.f, mwis_nx.f)


class TestMWCAll(unittest.TestCase):
    def test_three_approaches(self):
        for density in [np.random.random() for _ in range(5)]:
            for num_vertices in [np.random.randint(0, 30) for _ in range(5)]:
                for seed in range(5):
                    with self.subTest(
                        num_vertices=num_vertices, density=density, seed=seed
                    ):
                        num_vertices, density, seed = 10, density, seed
                        graph_sp = get_graph(num_vertices, density, seed, "spmatrix")
                        graph_pd = get_graph(num_vertices, density, seed, "DataFrame")
                        weights = np.random.randint(1, 100, size=num_vertices)
                        weights_pd = pd.DataFrame(weights, columns=["weights"])

                        mwc_sp = maximum_weighted_clique(graph_sp, weights)
                        mwc_pd = maximum_weighted_clique(graph_pd, weights_pd)
                        self.assertEqual(mwc_sp.f, mwc_pd.f)

                        if nx is not None:
                            graph_nx = get_graph(num_vertices, density, seed, "Graph")
                            mwc_nx = maximum_weighted_clique(graph_nx, weights)
                            self.assertEqual(mwc_nx.f, mwc_sp.f)
