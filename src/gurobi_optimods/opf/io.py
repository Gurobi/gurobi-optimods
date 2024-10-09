"""
Contains file I/O methods for reading and writing MATPOWER format.

- read_case_matpower: read case from MATPOWER format .mat file
- write_case_matpower: write case (or solution) to MATPOWER format .mat file
"""

import numpy as np
import pandas as pd
import scipy
import scipy.io

bus_field_names = [
    "bus_i",
    "type",
    "Pd",
    "Qd",
    "Gs",
    "Bs",
    "area",
    "Vm",
    "Va",
    "baseKV",
    "zone",
    "Vmax",
    "Vmin",
]
gen_field_names = [
    "bus",
    "Pg",
    "Qg",
    "Qmax",
    "Qmin",
    "Vg",
    "mBase",
    "status",
    "Pmax",
    "Pmin",
    "Pc1",
    "Pc2",
    "Qc1min",
    "Qc1max",
    "Qc2min",
    "Qc2max",
    "ramp_agc",
    "ramp_10",
    "ramp_30",
    "ramp_q",
    "apf",
]
branch_field_names = [
    "fbus",
    "tbus",
    "r",
    "x",
    "b",
    "rateA",
    "rateB",
    "rateC",
    "ratio",
    "angle",
    "status",
    "angmin",
    "angmax",
]
gencost_field_names = ["costtype", "startup", "shutdown", "n", "costvector"]


def read_case_matpower(file_path):
    """Read in a MATPOWER case file.

    The file must contain a MATPOWER version 2 'mpc' struct.

    Parameters
    ----------
    file_path : dict
        Path to .mat file

    Returns
    -------
    dict
        Case dictionary following MATPOWER notation
    """
    mat = scipy.io.loadmat(
        file_path, variable_names=["mpc"], struct_as_record=True, squeeze_me=True
    )

    if "mpc" not in mat.keys():
        raise ValueError("Provided .mat file does not have an mpc field")

    mpc = mat["mpc"]

    if mpc is None or int(mpc["version"]) != 2:
        raise ValueError("Provided .mat file must use MATPOWER specification version 2")

    missing_keys = sorted(
        {"baseMVA", "bus", "gen", "branch", "gencost"} - set(mpc.dtype.names)
    )
    if missing_keys:
        raise ValueError(f"Provided .mat file is missing keys: {missing_keys}")

    def fix_shape(arr, ncols=None):
        if ncols is None:
            ncols = arr.shape[-1]
        if arr.ndim == 1:
            return arr[:ncols].reshape((1, -1))
        else:
            assert arr.ndim == 2
            return arr[:, :ncols]

    # We assume column ordering matches the field name order, and discard any
    # additional columns.
    bus_array = fix_shape(mpc["bus"].item(), ncols=len(bus_field_names))
    bus = pd.DataFrame(data=bus_array, columns=bus_field_names).assign(
        bus_i=lambda df: df["bus_i"].astype(int),
        type=lambda df: df["type"].astype(int),
    )
    gen_array = fix_shape(mpc["gen"].item(), ncols=len(gen_field_names))

    if gen_array.shape[1] < len(gen_field_names):
        # If gen_array doesn't have all columns fill the remainder with Nan
        # We must have at least 10 columns
        assert gen_array.shape[1] >= 10
        missing = len(gen_field_names) - gen_array.shape[1]
        nan_array = np.full((gen_array.shape[0], missing), np.nan)
        gen_array = np.append(gen_array, nan_array, axis=1)
    gen = pd.DataFrame(data=gen_array, columns=gen_field_names).assign(
        bus=lambda df: df["bus"].astype(int)
    )

    branch_array = fix_shape(mpc["branch"].item(), ncols=len(branch_field_names))
    branch = pd.DataFrame(data=branch_array, columns=branch_field_names).assign(
        fbus=lambda df: df["fbus"].astype(int), tbus=lambda df: df["tbus"].astype(int)
    )

    # gencost has variable number of columns, we contract the extra columns to a
    # sub-list for costvector data. The length of this list should match the 'n'
    # field. Note: we don't check this here, it will be checked when the model
    # is built.
    gencost_array = fix_shape(mpc["gencost"].item())
    gencost_main = pd.DataFrame(
        data=gencost_array[:, :4], columns=gencost_field_names[:4]
    )
    gencost_vectors = gencost_array[:, 4:].tolist()
    gencost = [
        dict(main, costvector=vector)
        for main, vector in zip(gencost_main.to_dict("records"), gencost_vectors)
    ]

    return {
        "baseMVA": float(mpc["baseMVA"]),
        "bus": bus.to_dict("records"),
        "gen": gen.to_dict("records"),
        "branch": branch.to_dict("records"),
        "gencost": gencost,
        "casename": mpc["casename"] if "casename" in mpc.dtype.names else "case",
    }


def write_case_matpower(case, file_path):
    """Write a .mat file containing a single struct (named 'mpc') in MATPOWER
    format. 'case' follows the dictionary-based format we use for the solver API.
    """
    bus = pd.DataFrame(case["bus"])
    gen = pd.DataFrame(case["gen"])
    branch = pd.DataFrame(case["branch"])
    gencost = pd.DataFrame(case["gencost"])

    assert set(bus.columns) == set(bus_field_names)
    assert set(gen.columns) == set(gen_field_names)
    assert set(branch.columns) == set(branch_field_names)
    assert set(gencost.columns) == set(gencost_field_names)

    # Expand costvector list of lists to array
    part1 = gencost[gencost_field_names[:-1]].values
    part2 = np.array(list(gencost["costvector"].values))
    gencost_values = np.hstack([part1, part2])

    scipy.io.savemat(
        file_path,
        {
            "mpc": {
                "version": 2,
                "baseMVA": case["baseMVA"],
                "bus": bus[bus_field_names].values,
                "gen": gen[gen_field_names].values,
                "branch": branch[branch_field_names].values,
                "gencost": gencost_values,
                "casename": case.get("casename", "Undefined"),
            }
        },
    )
