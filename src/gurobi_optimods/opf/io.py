import numpy as np
import pandas as pd
import scipy

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


def read_case_matfile(file_path):
    """Read a .mat file containing a MATPOWER version 2 'mpc' struct. Returns
    the dictionary-based format we use for the solver APIs."""
    mat = scipy.io.loadmat(
        file_path, variable_names=["mpc"], struct_as_record=True, squeeze_me=True
    )
    mpc = mat["mpc"]
    assert int(mpc["version"]) == 2
    bus = pd.DataFrame(data=mpc["bus"].item(), columns=bus_field_names).assign(
        bus_i=lambda df: df["bus_i"].astype(int),
        type=lambda df: df["type"].astype(int),
    )
    gen = pd.DataFrame(data=mpc["gen"].item(), columns=gen_field_names).assign(
        bus=lambda df: df["bus"].astype(int)
    )
    branch = pd.DataFrame(data=mpc["branch"].item(), columns=branch_field_names).assign(
        fbus=lambda df: df["fbus"].astype(int), tbus=lambda df: df["tbus"].astype(int)
    )

    # gencost has variable number of columns, we contract the extra columns
    # to a sub-list for costvector data
    gencost_main = pd.DataFrame(
        data=mpc["gencost"].item()[:, :4], columns=gencost_field_names[:4]
    )
    gencost_vectors = mpc["gencost"].item()[:, 4:].tolist()
    gencost = [
        dict(main, costvector=vector)
        for main, vector in zip(gencost_main.to_dict("records"), gencost_vectors)
    ]
    for entry in gencost:
        assert entry["n"] == len(entry["costvector"])

    return {
        "baseMVA": float(mpc["baseMVA"]),
        "bus": bus.to_dict("records"),
        "gen": gen.to_dict("records"),
        "branch": branch.to_dict("records"),
        "gencost": gencost,
    }


def write_case_matfile(case, file_path):
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
            }
        },
    )
