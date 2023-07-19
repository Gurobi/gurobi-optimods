# Tests of reads/writes to data files

import json
import pathlib
import tempfile
import unittest

from gurobi_optimods.opf.io import read_case_matfile, write_case_matfile


class TestIO(unittest.TestCase):
    def test_roundtrip(self):
        # round trip test of dictionary format

        original = json.loads(
            pathlib.Path(__file__)
            .parent.joinpath("data/bus2-branch2-gen2.json")
            .read_text()
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpfile = pathlib.Path(tmpdir) / "testcase.mat"
            write_case_matfile(original, tmpfile)
            reread = read_case_matfile(tmpfile)

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

    def test_read_case(self):
        # Check that all example cases are read without errors

        from gurobi_optimods.datasets import load_caseopfmat

        case_mat_files = [
            load_caseopfmat(case) for case in ["9", "14", "57", "118", "300", "NY"]
        ]

        for file_path in case_mat_files:
            with self.subTest(file_path=file_path):
                # Should read without errors
                original = read_case_matfile(file_path)

                # Test write and read back
                with tempfile.TemporaryDirectory() as tmpdir:
                    tmpfile = pathlib.Path(tmpdir) / "testcase.mat"
                    write_case_matfile(original, tmpfile)
                    reread = read_case_matfile(tmpfile)

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
