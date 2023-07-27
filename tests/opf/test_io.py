# Tests of reads/writes to data files

import json
import pathlib
import tempfile
import unittest
from contextlib import contextmanager

import scipy

from gurobi_optimods.datasets import load_opf_example, load_opf_extra
from gurobi_optimods.opf import read_case_matpower, write_case_matpower


class TestBadData(unittest.TestCase):
    @contextmanager
    def save_scipy_data_to_tempfile(self, content):
        with tempfile.TemporaryDirectory() as tmpdir:
            matfile = pathlib.Path(tmpdir) / "testcase.mat"
            scipy.io.savemat(matfile, content)
            yield matfile

    def test_no_mpc_key(self):
        content = {}
        with self.save_scipy_data_to_tempfile(content) as matfile:
            with self.assertRaisesRegex(
                ValueError, "Provided .mat file does not have an mpc field"
            ):
                read_case_matpower(matfile)

    def test_bad_version(self):
        content = {"mpc": {}}
        with self.save_scipy_data_to_tempfile(content) as matfile:
            with self.assertRaisesRegex(
                ValueError,
                "Provided .mat file must use MATPOWER specification version 2",
            ):
                read_case_matpower(matfile)

        content = {"mpc": {"version": 1}}
        with self.save_scipy_data_to_tempfile(content) as matfile:
            with self.assertRaisesRegex(
                ValueError,
                "Provided .mat file must use MATPOWER specification version 2",
            ):
                read_case_matpower(matfile)

    def test_missing_field(self):
        content = {"mpc": {"version": 2, "branch": 3, "gen": 4, "baseMVA": 5}}
        with self.save_scipy_data_to_tempfile(content) as matfile:
            with self.assertRaisesRegex(
                ValueError, "Provided .mat file is missing keys"
            ):
                read_case_matpower(matfile)


class TestDatasets(unittest.TestCase):
    def test_load_opf_example(self):
        # Check that example cases can be loaded by name
        for case in ["case9", "case14", "case57", "case118", "case300", "caseNY"]:
            case_data = load_opf_example(case)
            self.assertEqual(
                set(case_data.keys()), {"baseMVA", "bus", "branch", "gen", "gencost"}
            )

    def test_load_opf_extras_coordinates(self):
        for extra in ["case9-coordinates", "caseNY-coordinates"]:
            data = load_opf_extra(extra)
            for bus, (lat, lon) in data.items():
                self.assertIsInstance(bus, int)
                self.assertIsInstance(lat, float)
                self.assertIsInstance(lon, float)

    def test_load_opf_extras_case9voltages(self):
        data = load_opf_extra("case9-voltages")
        self.assertEqual(data[2], (1.099999, 20.552543))
        for bus, (Vm, Va) in data.items():
            self.assertIsInstance(bus, int)
            self.assertIsInstance(Vm, float)
            self.assertIsInstance(Va, float)


class TestIO(unittest.TestCase):
    def setUp(self):
        here = pathlib.Path(__file__).parent.resolve()
        self.testdata_dir = here.joinpath("data")
        self.dataset_dir = here.parent.parent.joinpath("src/gurobi_optimods/data/opf")

    def test_matfile_roundtrip(self):
        # round trip test of dictionary format

        for jsondata in ["bus1-branch1-gen1.json", "bus2-branch2-gen2.json"]:
            with self.subTest(jsondata=jsondata):
                original = json.loads(self.testdata_dir.joinpath(jsondata).read_text())

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
