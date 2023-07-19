# Tests of reads/writes to data files

import json
import pathlib
import tempfile
import unittest

from gurobi_optimods.datasets import load_opf_example
from gurobi_optimods.opf import read_case_matpower, write_case_matpower


class TestDatasets(unittest.TestCase):
    def test_load_opf_example(self):
        # Check that example cases can be loaded by name
        for case in ["case9", "case14", "case57", "case118", "case300", "caseNY"]:
            case_data = load_opf_example(case)
            self.assertEqual(
                set(case_data.keys()), {"baseMVA", "bus", "branch", "gen", "gencost"}
            )


class TestIO(unittest.TestCase):
    def setUp(self):
        here = pathlib.Path(__file__).parent.resolve()
        self.testdata_dir = here.joinpath("data")
        self.dataset_dir = here.parent.parent.joinpath("src/gurobi_optimods/data/opf")

    def test_matfile_roundtrip(self):
        # round trip test of dictionary format

        original = json.loads(
            self.testdata_dir.joinpath("bus2-branch2-gen2.json").read_text()
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpfile = pathlib.Path(tmpdir) / "testcase.mat"
            write_case_matpower(original, tmpfile)
            reread = read_case_matpower(tmpfile)

        self.assertEqual(set(reread.keys()), set(original.keys()))

        for field, data in reread.items():
            self.assertEqual(data, original[field])

        # Check types for some special cases
        for bus in reread["bus"]:
            self.assertIsInstance(bus["bus_i"], int)
            self.assertIsInstance(bus["type"], int)
        for gen in reread["gen"]:
            self.assertIsInstance(gen["bus"], int)
        for branch in reread["branch"]:
            self.assertIsInstance(branch["fbus"], int)
            self.assertIsInstance(branch["tbus"], int)

    def test_read_case_matpower(self):
        # Check that all example cases are read without errors

        case_mat_files = [
            self.dataset_dir.joinpath(f"case{case}.mat")
            for case in ["9", "14", "57", "118", "300", "NY"]
        ]

        for file_path in case_mat_files:
            with self.subTest(file_path=file_path):
                # Should read without errors
                original = read_case_matpower(file_path)

                # Test write and read back
                with tempfile.TemporaryDirectory() as tmpdir:
                    tmpfile = pathlib.Path(tmpdir) / "testcase.mat"
                    write_case_matpower(original, tmpfile)
                    reread = read_case_matpower(tmpfile)

                # The first read and the round-trip should match exactly
                self.assertEqual(set(reread.keys()), set(original.keys()))
                for field, data in reread.items():
                    self.assertEqual(data, original[field])

                # Check types for some special cases
                for bus in reread["bus"]:
                    self.assertIsInstance(bus["bus_i"], int)
                    self.assertIsInstance(bus["type"], int)
                for gen in reread["gen"]:
                    self.assertIsInstance(gen["bus"], int)
                for branch in reread["branch"]:
                    self.assertIsInstance(branch["fbus"], int)
                    self.assertIsInstance(branch["tbus"], int)
