import unittest

import numpy as np
import pandas as pd
import scipy.sparse as sp

try:
    import networkx as nx
except ImportError:
    nx = None

from gurobi_optimods import datasets
from gurobi_optimods.hamiltonian_cycle import min_bottleneck_hamiltonian_cycle


def frame(arcs):
    """Build the pandas input from {(source, target): weight}."""
    return pd.DataFrame(
        [(s, t, w) for (s, t), w in arcs.items()],
        columns=["source", "target", "weight"],
    ).set_index(["source", "target"])


def as_at(season, played_through):
    """The season as it stood at the end of ``played_through``: later games are
    on the fixture, but their results are not known yet."""
    return season.assign(winner=season["winner"].mask(season["round"] > played_through))


def build_arcs(season, by="round"):
    """The "u has beaten v" arcs implied by a season. A game the season records
    no winner for has not been decided yet, so it could go either way and
    contributes both orientations. Arcs are weighted by the ``by`` column of
    the game that created them. This mirrors the transform shown in the
    documentation."""
    undecided = season["winner"].isna()
    decided = season[~undecided & (season["winner"] != "draw")]
    loser = decided["away"].where(decided["winner"] == decided["home"], decided["home"])
    known = pd.DataFrame(
        {"source": decided["winner"], "target": loser, "weight": decided[by]}
    )

    later = season[undecided]
    possible = pd.concat(
        [
            pd.DataFrame(
                {
                    "source": later["home"],
                    "target": later["away"],
                    "weight": later[by],
                }
            ),
            pd.DataFrame(
                {
                    "source": later["away"],
                    "target": later["home"],
                    "weight": later[by],
                }
            ),
        ]
    )

    return (
        pd.concat([known, possible])
        .groupby(["source", "target"])["weight"]
        .min()
        .to_frame()
    )


class CycleAssertions:
    """Solutions are rarely unique, so validate the structure of whatever
    cycle comes back rather than pinning one particular answer."""

    def assert_valid_cycle(self, result, arcs, num_nodes):
        self.assertEqual(len(result.cycle), num_nodes)
        self.assertEqual(len(set(result.cycle)), num_nodes)
        self.assertEqual(len(result.weights), num_nodes)

        for position, node in enumerate(result.cycle):
            successor = result.cycle[(position + 1) % num_nodes]
            self.assertIn(
                (node, successor), arcs, f"{node} -> {successor} is not an arc"
            )
            self.assertEqual(result.weights[position], arcs[node, successor])

        self.assertEqual(result.bottleneck, max(result.weights))


class TestMinBottleneckHamiltonianCycle(unittest.TestCase, CycleAssertions):
    def test_triangle(self):
        arcs = {("a", "b"): 1, ("b", "c"): 2, ("c", "a"): 3}
        result = min_bottleneck_hamiltonian_cycle(frame(arcs), verbose=False)
        self.assert_valid_cycle(result, arcs, 3)
        self.assertEqual(result.bottleneck, 3)

    def test_min_max_is_not_min_sum(self):
        # Two cycles exist. a->b->c->d->a has weights {1, 6, 2, 3}: max 6,
        # sum 12. a->c->b->d->a has weights {4, 3, 5, 3}: max 5, sum 15. The
        # bottleneck objective must take the second, a min-sum one the first.
        arcs = {
            ("a", "b"): 1,
            ("b", "c"): 6,
            ("c", "d"): 2,
            ("d", "a"): 3,
            ("a", "c"): 4,
            ("c", "b"): 3,
            ("b", "d"): 5,
        }
        result = min_bottleneck_hamiltonian_cycle(frame(arcs), verbose=False)
        self.assert_valid_cycle(result, arcs, 4)
        self.assertEqual(result.bottleneck, 5)
        self.assertEqual(sum(result.weights), 15)

    def test_no_hamiltonian_cycle(self):
        # Every node has an arc in and an arc out, but a->b is a bridge that
        # cannot be returned across, so no cycle covers all four nodes.
        arcs = {
            ("a", "b"): 1,
            ("b", "a"): 1,
            ("c", "d"): 1,
            ("d", "c"): 1,
            ("b", "c"): 1,
        }
        with self.assertRaisesRegex(ValueError, "no Hamiltonian cycle"):
            min_bottleneck_hamiltonian_cycle(frame(arcs), verbose=False)

    def test_node_without_outgoing_arc_is_named(self):
        arcs = {("a", "b"): 1, ("b", "c"): 2, ("c", "a"): 3, ("a", "sink"): 4}
        with self.assertRaisesRegex(ValueError, "'sink' has no outgoing arcs"):
            min_bottleneck_hamiltonian_cycle(frame(arcs), verbose=False)

    def test_node_without_incoming_arc_is_named(self):
        arcs = {("a", "b"): 1, ("b", "c"): 2, ("c", "a"): 3, ("source", "a"): 4}
        with self.assertRaisesRegex(ValueError, "'source' has no incoming arcs"):
            min_bottleneck_hamiltonian_cycle(frame(arcs), verbose=False)

    def test_numeric_node_is_named_plainly(self):
        # Numeric labels arrive from a dataframe as numpy scalars, but the
        # message should read as the user wrote them.
        arcs = {(0, 1): 1, (1, 2): 2, (2, 0): 3, (0, 9): 4}
        with self.assertRaisesRegex(ValueError, r"Node 9 has no outgoing arcs"):
            min_bottleneck_hamiltonian_cycle(frame(arcs), verbose=False)

    def test_too_few_nodes(self):
        arcs = {("a", "b"): 1, ("b", "a"): 2}
        with self.assertRaisesRegex(ValueError, "at least 3 nodes"):
            min_bottleneck_hamiltonian_cycle(frame(arcs), verbose=False)

    def test_self_loops_ignored(self):
        arcs = {("a", "b"): 1, ("b", "c"): 2, ("c", "a"): 3}
        with_loops = {**arcs, ("a", "a"): 0, ("b", "b"): 0}
        result = min_bottleneck_hamiltonian_cycle(frame(with_loops), verbose=False)
        self.assert_valid_cycle(result, arcs, 3)
        self.assertEqual(result.bottleneck, 3)

    def test_parallel_arcs_keep_earliest(self):
        # a->b appears twice; only the earlier weight should ever be used.
        rows = pd.DataFrame(
            [("a", "b", 9), ("a", "b", 1), ("b", "c", 2), ("c", "a", 3)],
            columns=["source", "target", "weight"],
        ).set_index(["source", "target"])
        result = min_bottleneck_hamiltonian_cycle(rows, verbose=False)
        self.assertEqual(result.bottleneck, 3)
        self.assertEqual(sorted(result.weights), [1, 2, 3])

    def test_weights_are_optional(self):
        rows = pd.DataFrame(
            [("a", "b"), ("b", "c"), ("c", "a")], columns=["source", "target"]
        ).set_index(["source", "target"])
        result = min_bottleneck_hamiltonian_cycle(rows, verbose=False)
        self.assert_valid_cycle(
            result, {k: 1 for k in [("a", "b"), ("b", "c"), ("c", "a")]}, 3
        )
        self.assertEqual(result.bottleneck, 1)

    def test_datetime_weights(self):
        arcs = {
            ("a", "b"): pd.Timestamp("2026-03-05"),
            ("b", "c"): pd.Timestamp("2026-04-11"),
            ("c", "a"): pd.Timestamp("2026-03-22"),
        }
        result = min_bottleneck_hamiltonian_cycle(frame(arcs), verbose=False)
        self.assert_valid_cycle(result, arcs, 3)
        self.assertEqual(result.bottleneck, pd.Timestamp("2026-04-11"))
        self.assertIsInstance(result.bottleneck, pd.Timestamp)

    def test_node_labels_are_plain_python(self):
        # As with the weights, labels come back in their natural python form
        # rather than as the numpy scalars a dataframe stores them in.
        arcs = {(0, 1): 1, (1, 2): 2, (2, 0): 3}
        result = min_bottleneck_hamiltonian_cycle(frame(arcs), verbose=False)
        self.assert_valid_cycle(result, arcs, 3)
        for node in result.cycle:
            self.assertIs(type(node), int)

    def test_source_target_columns_accepted(self):
        rows = pd.DataFrame(
            [("a", "b", 1), ("b", "c", 2), ("c", "a", 3)],
            columns=["source", "target", "weight"],
        )
        result = min_bottleneck_hamiltonian_cycle(rows, verbose=False)
        self.assertEqual(result.bottleneck, 3)

    def test_bad_frame_layout(self):
        rows = pd.DataFrame([("a", "b", 1)], columns=["from", "to", "weight"])
        with self.assertRaisesRegex(ValueError, "source"):
            min_bottleneck_hamiltonian_cycle(rows, verbose=False)

    def test_unknown_graph_type(self):
        with self.assertRaisesRegex(ValueError, "Unknown graph type"):
            min_bottleneck_hamiltonian_cycle([(0, 1), (1, 0)], verbose=False)


class TestInputTypes(unittest.TestCase, CycleAssertions):
    """The same graph, expressed three ways, must give the same answer."""

    def setUp(self):
        self.arcs = {
            (0, 1): 1,
            (1, 2): 6,
            (2, 3): 2,
            (3, 0): 3,
            (0, 2): 4,
            (2, 1): 3,
            (1, 3): 5,
        }

    def test_pandas(self):
        result = min_bottleneck_hamiltonian_cycle(frame(self.arcs), verbose=False)
        self.assert_valid_cycle(result, self.arcs, 4)
        self.assertEqual(result.bottleneck, 5)

    def test_scipy(self):
        rows = np.array([s for s, _ in self.arcs])
        cols = np.array([t for _, t in self.arcs])
        data = np.array(list(self.arcs.values()))
        matrix = sp.coo_array((data, (rows, cols)), shape=(4, 4))
        result = min_bottleneck_hamiltonian_cycle(matrix, verbose=False)
        self.assert_valid_cycle(result, self.arcs, 4)
        self.assertEqual(result.bottleneck, 5)

    def test_scipy_keeps_zero_weights(self):
        # An arc of weight zero is a stored entry like any other, and must not
        # be read as an absent arc.
        arcs = {(0, 1): 0, (1, 2): 2, (2, 0): 1}
        rows = np.array([s for s, _ in arcs])
        cols = np.array([t for _, t in arcs])
        data = np.array(list(arcs.values()))
        matrix = sp.coo_array((data, (rows, cols)), shape=(3, 3))
        result = min_bottleneck_hamiltonian_cycle(matrix, verbose=False)
        self.assert_valid_cycle(result, arcs, 3)
        self.assertEqual(result.bottleneck, 2)

    def test_scipy_non_square(self):
        matrix = sp.coo_array(([1, 2], ([0, 1], [1, 2])), shape=(3, 4))
        with self.assertRaisesRegex(ValueError, "square"):
            min_bottleneck_hamiltonian_cycle(matrix, verbose=False)

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_networkx(self):
        graph = nx.DiGraph()
        for (source, target), weight in self.arcs.items():
            graph.add_edge(source, target, weight=weight)
        result = min_bottleneck_hamiltonian_cycle(graph, verbose=False)
        self.assert_valid_cycle(result, self.arcs, 4)
        self.assertEqual(result.bottleneck, 5)

    @unittest.skipIf(nx is None, "networkx is not installed")
    def test_networkx_undirected_traversable_both_ways(self):
        # An undirected edge may be used in either direction, so the square
        # graph has a Hamiltonian cycle even though no orientation is given.
        graph = nx.Graph()
        graph.add_edge("a", "b", weight=1)
        graph.add_edge("b", "c", weight=2)
        graph.add_edge("c", "d", weight=3)
        graph.add_edge("d", "a", weight=4)
        result = min_bottleneck_hamiltonian_cycle(graph, verbose=False)
        self.assertEqual(len(result.cycle), 4)
        self.assertEqual(result.bottleneck, 4)


class TestAFLSeason(unittest.TestCase, CycleAssertions):
    def setUp(self):
        self.season = datasets.load_afl_season()

    def test_dataset(self):
        self.assertEqual(len(self.season), 207)
        self.assertEqual((self.season["winner"] == "draw").sum(), 3)
        teams = set(self.season["home"]) | set(self.season["away"])
        self.assertEqual(len(teams), 18)
        self.assertEqual(self.season["round"].min(), 0)
        self.assertEqual(self.season["round"].max(), 24)
        self.assertEqual(self.season["date"].max(), pd.Timestamp("2026-08-23"))

    def test_completed_season(self):
        arcs = build_arcs(self.season)
        result = min_bottleneck_hamiltonian_cycle(arcs, verbose=False)
        self.assert_valid_cycle(result, arcs["weight"].to_dict(), 18)
        # No team can appear in the cycle before it has both a win and a loss,
        # and RIC did not win until round 8.
        self.assertGreaterEqual(result.bottleneck, 8)
        self.assertEqual(result.bottleneck, 8)

    def test_preseason_fixture_only(self):
        # Nothing played: every game offers both orientations.
        arcs = build_arcs(as_at(self.season, -1))
        result = min_bottleneck_hamiltonian_cycle(arcs, verbose=False)
        self.assert_valid_cycle(result, arcs["weight"].to_dict(), 18)
        # Each team plays once per round, so two rounds are needed at minimum.
        self.assertGreaterEqual(result.bottleneck, 2)
        self.assertEqual(result.bottleneck, 3)

    def test_answer_is_monotone_in_results_known(self):
        # Learning a result only ever removes arcs, so the earliest possible
        # completion can never move earlier as the season progresses.
        answers = [
            min_bottleneck_hamiltonian_cycle(
                build_arcs(as_at(self.season, k)), verbose=False
            ).bottleneck
            for k in range(-1, 10)
        ]
        self.assertEqual(answers, sorted(answers))
        self.assertEqual(answers[0], 3)
        self.assertEqual(answers[-1], 8)

    def test_season_in_progress(self):
        # A season still being played records no winner for its later games,
        # and every one of those contributes both orientations.
        arcs = build_arcs(as_at(self.season, 5))
        undecided = self.season[self.season["round"] > 5]
        for home, away in zip(undecided["home"], undecided["away"]):
            self.assertIn((home, away), arcs.index)
            self.assertIn((away, home), arcs.index)
        self.assertGreater(len(arcs), len(build_arcs(self.season)))

    def test_weighting_by_date(self):
        arcs = build_arcs(self.season, by="date")
        result = min_bottleneck_hamiltonian_cycle(arcs, verbose=False)
        self.assert_valid_cycle(result, arcs["weight"].to_dict(), 18)
        self.assertIsInstance(result.bottleneck, pd.Timestamp)
        # Dates and rounds order the games the same way, so the answer must
        # land inside the round the round-weighted model picked.
        rounds = self.season.loc[self.season["round"] == 8, "date"]
        self.assertTrue(rounds.min() <= result.bottleneck <= rounds.max())


if __name__ == "__main__":
    unittest.main()
