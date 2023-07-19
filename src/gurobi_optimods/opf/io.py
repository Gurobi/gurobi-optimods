"""
Contains file I/O methods for reading and writing known data formats.

- read_case_matfile: read case from MATPOWER format .mat file
- write_case_matfile: write case (or solution) to MATPOWER format .mat file
- read_coords_from_csv_file: read coordinate data
- read_voltages_from_csv_file: read voltage data
"""

import csv

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


"""
Some errors thrown by the old reader (TODO add tests):

    if "mpc" not in mat.keys():
        raise ValueError("Provided .mat file does not have an mpc field")

    for x in ["baseMVA", "bus", "gen", "branch", "gencost"]:
        if x not in mpc.keys():
            raise ValueError(f"Provided .mat file does not have a {x} field")

    # TODO this case should be handled? Also still need to handle 1gen, 1branch,
    etc if mpcbuses.ndim == 1:
        raise ValueError("Provided .mat files has only 1 bus")

    if numgencosts > numgens:
        # FIXME: spec says we can have twice as many entries, representing
        reactive power. Do we not handle this case?

        raise ValueError(f"Read {numgencosts} gen costs but only {numgens}
        generators.")
"""


def read_case_matfile(file_path):
    """Read a .mat file containing a MATPOWER version 2 'mpc' struct. Returns
    the dictionary-based format we use for the solver APIs."""
    mat = scipy.io.loadmat(
        file_path, variable_names=["mpc"], struct_as_record=True, squeeze_me=True
    )
    mpc = mat["mpc"]
    assert int(mpc["version"]) == 2

    # We assume column ordering matches the field name order, and discard any
    # additional columns.
    bus_array = mpc["bus"].item()
    bus_array = bus_array[:, : len(bus_field_names)]
    bus = pd.DataFrame(data=bus_array, columns=bus_field_names).assign(
        bus_i=lambda df: df["bus_i"].astype(int),
        type=lambda df: df["type"].astype(int),
    )
    gen_array = mpc["gen"].item()
    gen_array = gen_array[:, : len(gen_field_names)]
    gen = pd.DataFrame(data=gen_array, columns=gen_field_names).assign(
        bus=lambda df: df["bus"].astype(int)
    )
    branch_array = mpc["branch"].item()
    branch_array = branch_array[:, : len(branch_field_names)]
    branch = pd.DataFrame(data=branch_array, columns=branch_field_names).assign(
        fbus=lambda df: df["fbus"].astype(int), tbus=lambda df: df["tbus"].astype(int)
    )

    # gencost has variable number of columns, we contract the extra columns
    # to a sub-list for costvector data. The length of this list should match
    # the 'n' field.
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


def read_file_csv(filename, data):
    """
    Reads bus data from a `.csv` file

    :param filename: Path to text based file holding bus data
    :type filename: str
    :param data: Name of data to be read in
    :type data: str

    :return: Dictionary holding the respective data of the form busID : (data1, data2)
    :rtype: dict
    """

    with open(filename, mode="r") as infile:
        reader = csv.reader(infile)
        # Sanity check
        data_dict = {}
        first = True
        for rows in reader:
            # Skip first row of the .csv file
            if first:
                first = False
                continue

            if len(rows) != 5:
                # TODO ... what???
                raise ValueError(
                    f"Incorrect input in {data} .csv file {filename}. Number of columns does not equal 5."
                )

            data_dict[int(rows[1])] = (
                float(rows[3]),
                float(rows[4]),
            )

    return data_dict


def turn_opf_dict_into_mat_file(solution, filename):
    """
    Writes a `.mat` file out of an OPF solution dictionary

    :param solution: OPF solution dictionary
    :type solution: dict
    :param filename: Name of `.mat` file where to write the solution data
    :type filename: str

    FIXME this is not used anywhere, and cannot possibly work. It just turns
    a dict of dicts into an array with no care for the keys (the order cannot
    be assumed in general).

    write_case_matfile can do this, but there are additional fields we add not
    specified in matpower, so its not clear what order they should be written in
    """

    # Buses
    buses = solution["bus"]
    matrix = []
    for bus in buses.values():
        matrix.append(list(bus.values()))

    solution["bus"] = np.array(matrix)
    # Generators
    gens = solution["gen"]
    matrix = []
    for gen in gens.values():
        matrix.append(list(gen.values()))

    solution["gen"] = np.array(matrix)
    # Branches
    branches = solution["branch"]
    matrix = []
    for branch in branches.values():
        matrix.append(list(branch.values()))

    solution["branch"] = np.array(matrix)
    # Generator costs
    gencosts = solution["gencost"]
    matrix = []
    # For generator costs, the last dictionary entry is a list of values
    # TODO NOT SAFE TO RELY ON ORDER
    for genc in gencosts.values():
        l = list(genc.values())
        costvector = l[-1]
        l.pop()
        for item in costvector:
            l.append(item)

        matrix.append(l)

    solution["gencost"] = np.array(matrix)
    # Write mat file
    scipy.io.savemat(filename, {"result": solution})


def read_coords_from_csv_file(coordsfile):
    """
    Helper function for users. Constructs a coordinate dictionary which can be
    used as input for the ``generate_opf_solution_figure`` function.

    :param coordsfile: Name of and possibly full path to bus coordinates file
        given as `.csv` file. The `.csv` file has to consist of the following
        columns in the given order: ``index(starting with 0), busID, busname,
        latitude, longitude``
    :type coordsfile: str

    :return: Dictionary of the given coordinates which can be used in the
             ``generate_opf_solution_figure`` function
    :rtype: dict

    .. note::
        The first row of the csv file is ignored.
    """

    coord_dict = read_file_csv(coordsfile, "coordinates")

    return coord_dict


def read_voltages_from_csv_file(voltsfile):
    """
    Helper function for users. Constructs a bus input voltage dictionary which
    can be used as input for the ``check_voltage_solution_violations`` function

    :param voltsfile: Name of and possibly full path to voltage input file given
        as `.csv` file. The `.csv` file has to consist of the following columns
        in the given order: ``index(starting with 0), busID, busname, latitude,
        longitude``
    :type voltsfile: str

    :return: Dictionary of the given coordinates which can be used in the
             ``check_voltage_solution_violations`` function
    :rtype: dict

    .. note::
        The first row of the csv file is ignored.
    """

    volts_dict = read_file_csv(voltsfile, "voltages")

    return volts_dict
