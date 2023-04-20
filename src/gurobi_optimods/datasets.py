"""
Module for loading datasets for use in optimods examples, in the same vein
as sklearn.datasets.
"""

import pathlib

import pandas as pd

from collections import OrderedDict

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
            DATA_FILE_DIR / "workforce/availability.csv", parse_dates=["Shift"]
        ),
        pay_rates=pd.read_csv(DATA_FILE_DIR / "workforce/pay_rates.csv"),
        shift_requirements=pd.read_csv(
            DATA_FILE_DIR / "workforce/shift_requirements.csv", parse_dates=["Shift"]
        ),
    )


def load_opfsettings():
    conf = str(DATA_FILE_DIR) + "/opf/opfsettings.txt"
    return conf


def load_caseopf(number):
    case = str(DATA_FILE_DIR) + "/opf/case" + number + ".m"
    return case


def load_coordsfilepath(filename):
    file = str(DATA_FILE_DIR) + "/opf/" + filename
    return file


def load_caseopfmat(number):
    case = str(DATA_FILE_DIR) + "/opf/case" + number + ".mat"
    return case


def load_caseNYopf():  # real world data case
    case = str(DATA_FILE_DIR) + "/opf/caseNY.m"
    casemat = str(DATA_FILE_DIR) + "/opf/caseNY.mat"
    return case, casemat


def load_opfgraphicssettings():
    conf = str(DATA_FILE_DIR) + "/opf/graphicssettings.txt"
    return conf


def load_opfdictcase():
    casefile_dict = {
        "baseMVA": 100.0,
        "bus": {
            1: {
                "count": 1,  # TODO-Dan Do we need the count? If yes, why?
                "bus_i": 1,
                "type": 3,
                "Pd": 0.0,
                "Qd": 0.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "area": 0.0,
                "Vm": 1.0,
                "Va": 1.0,
                "baseKV": 345.0,
                "zone": 1.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
                "lnum": 18,  # TODO-Dan Do we need the line number? If yes, why?
            },
            2: {
                "count": 2,
                "bus_i": 2,
                "type": 2,
                "Pd": 0.0,
                "Qd": 0.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "area": 0.0,
                "Vm": 1.0,
                "Va": 1.0,
                "baseKV": 345.0,
                "zone": 1.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
                "lnum": 19,
            },
            3: {
                "count": 3,
                "bus_i": 3,
                "type": 2,
                "Pd": 0.0,
                "Qd": 0.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "area": 0.0,
                "Vm": 1.0,
                "Va": 1.0,
                "baseKV": 345.0,
                "zone": 1.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
                "lnum": 20,
            },
            4: {
                "count": 4,
                "bus_i": 4,
                "type": 1,
                "Pd": 0.0,
                "Qd": 0.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "area": 0.0,
                "Vm": 1.0,
                "Va": 1.0,
                "baseKV": 345.0,
                "zone": 1.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
                "lnum": 21,
            },
            5: {
                "count": 5,
                "bus_i": 5,
                "type": 1,
                "Pd": 90.0,
                "Qd": 30.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "area": 0.0,
                "Vm": 1.0,
                "Va": 1.0,
                "baseKV": 345.0,
                "zone": 1.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
                "lnum": 22,
            },
            6: {
                "count": 6,
                "bus_i": 6,
                "type": 1,
                "Pd": 0.0,
                "Qd": 0.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "area": 0.0,
                "Vm": 1.0,
                "Va": 1.0,
                "baseKV": 345.0,
                "zone": 1.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
                "lnum": 23,
            },
            7: {
                "count": 7,
                "bus_i": 7,
                "type": 1,
                "Pd": 100.0,
                "Qd": 35.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "area": 0.0,
                "Vm": 1.0,
                "Va": 1.0,
                "baseKV": 345.0,
                "zone": 1.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
                "lnum": 24,
            },
            8: {
                "count": 8,
                "bus_i": 8,
                "type": 1,
                "Pd": 0.0,
                "Qd": 0.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "area": 0.0,
                "Vm": 1.0,
                "Va": 1.0,
                "baseKV": 345.0,
                "zone": 1.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
                "lnum": 25,
            },
            9: {
                "count": 9,
                "bus_i": 9,
                "type": 1,
                "Pd": 125.0,
                "Qd": 50.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "area": 0.0,
                "Vm": 1.0,
                "Va": 1.0,
                "baseKV": 345.0,
                "zone": 1.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
                "lnum": 26,
            },
        },
        "gen": {
            1: {
                "gencount1": 1,  # TODO-Dan Do we need the gencount1 value? If yes, why?
                "bus": 1,
                "Pg": 0.0,
                "Qg": 0.0,
                "Qmax": 300.0,
                "Qmin": -300.0,
                "Vg": 1,
                "mBase": 100,
                "status": 1,
                "Pmax": 250.0,
                "Pmin": 10.0,
                "Pc1": 0,
                "Pc2": 0,
                "Qc1min": 0,
                "Qc1max": 0,
                "Qc2min": 0,
                "Qc2max": 0,
                "ramp_agc": 0,
                "ramp_10": 0,
                "ramp_30": 0,
                "ramp_q": 0,
                "apf": 0,
                "lnum": 32,  # TODO-Dan Do we need the line number? If yes, why?
            },
            2: {
                "gencount1": 2,
                "bus": 2,
                "Pg": 163.0,
                "Qg": 0.0,
                "Qmax": 300.0,
                "Qmin": -300.0,
                "Vg": 1,
                "mBase": 100,
                "status": 1,
                "Pmax": 300.0,
                "Pmin": 10.0,
                "Pc1": 0,
                "Pc2": 0,
                "Qc1min": 0,
                "Qc1max": 0,
                "Qc2min": 0,
                "Qc2max": 0,
                "ramp_agc": 0,
                "ramp_10": 0,
                "ramp_30": 0,
                "ramp_q": 0,
                "apf": 0,
                "lnum": 33,
            },
            3: {
                "gencount1": 3,
                "bus": 3,
                "Pg": 85.0,
                "Qg": 0.0,
                "Qmax": 300.0,
                "Qmin": -300.0,
                "Vg": 1,
                "mBase": 100,
                "status": 1,
                "Pmax": 270.0,
                "Pmin": 10.0,
                "Pc1": 0,
                "Pc2": 0,
                "Qc1min": 0,
                "Qc1max": 0,
                "Qc2min": 0,
                "Qc2max": 0,
                "ramp_agc": 0,
                "ramp_10": 0,
                "ramp_30": 0,
                "ramp_q": 0,
                "apf": 0,
                "lnum": 34,
            },
        },
        "branch": {
            1: {
                "branchcount1": 1,  # TODO-Dan Do we need the branchcount? If yes, why?
                "fbus": 1,
                "tbus": 4,
                "r": 0.0,
                "x": 0.0576,
                "b": 0.0,
                "rateA": 250.0,
                "rateB": 250.0,
                "rateC": 250.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "angmin": -360.0,
                "angmax": 360.0,
                "lnum": 40,  # TODO-Dan Do we need the line number? If yes, why?
            },
            2: {
                "branchcount1": 2,
                "fbus": 4,
                "tbus": 5,
                "r": 0.017,
                "x": 0.092,
                "b": 0.158,
                "rateA": 250.0,
                "rateB": 250.0,
                "rateC": 250.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "angmin": -360.0,
                "angmax": 360.0,
                "lnum": 41,
            },
            3: {
                "branchcount1": 3,
                "fbus": 5,
                "tbus": 6,
                "r": 0.039,
                "x": 0.17,
                "b": 0.358,
                "rateA": 150.0,
                "rateB": 150.0,
                "rateC": 150.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "angmin": -360.0,
                "angmax": 360.0,
                "lnum": 42,
            },
            4: {
                "branchcount1": 4,
                "fbus": 3,
                "tbus": 6,
                "r": 0.0,
                "x": 0.0586,
                "b": 0.0,
                "rateA": 300.0,
                "rateB": 300.0,
                "rateC": 300.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "angmin": -360.0,
                "angmax": 360.0,
                "lnum": 43,
            },
            5: {
                "branchcount1": 5,
                "fbus": 6,
                "tbus": 7,
                "r": 0.0119,
                "x": 0.1008,
                "b": 0.209,
                "rateA": 150.0,
                "rateB": 150.0,
                "rateC": 150.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "angmin": -360.0,
                "angmax": 360.0,
                "lnum": 44,
            },
            6: {
                "branchcount1": 6,
                "fbus": 7,
                "tbus": 8,
                "r": 0.0085,
                "x": 0.072,
                "b": 0.149,
                "rateA": 250.0,
                "rateB": 250.0,
                "rateC": 250.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "angmin": -360.0,
                "angmax": 360.0,
                "lnum": 45,
            },
            7: {
                "branchcount1": 7,
                "fbus": 8,
                "tbus": 2,
                "r": 0.0,
                "x": 0.0625,
                "b": 0.0,
                "rateA": 250.0,
                "rateB": 250.0,
                "rateC": 250.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "angmin": -360.0,
                "angmax": 360.0,
                "lnum": 46,
            },
            8: {
                "branchcount1": 8,
                "fbus": 8,
                "tbus": 9,
                "r": 0.032,
                "x": 0.161,
                "b": 0.306,
                "rateA": 250.0,
                "rateB": 250.0,
                "rateC": 250.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "angmin": -360.0,
                "angmax": 360.0,
                "lnum": 47,
            },
            9: {
                "branchcount1": 9,
                "fbus": 9,
                "tbus": 4,
                "r": 0.01,
                "x": 0.085,
                "b": 0.176,
                "rateA": 250.0,
                "rateB": 250.0,
                "rateC": 250.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "angmin": -360.0,
                "angmax": 360.0,
                "lnum": 48,
            },
        },
        "gencost": {
            1: {
                "costtype": 2,
                "startup": 1500,
                "shutdown": 0,
                "n": 3,
                "costvector": [0.11, 5, 150],
            },
            2: {
                "costtype": 2,
                "startup": 2000,
                "shutdown": 0,
                "n": 3,
                "costvector": [0.085, 1.2, 600],
            },
            3: {
                "costtype": 2,
                "startup": 3000,
                "shutdown": 0,
                "n": 3,
                "costvector": [0.1225, 1, 335],
            },
        },
    }
    return casefile_dict


def load_case9solution():
    solution = OrderedDict(
        [
            ("theta_1", 6.177764273002596),
            ("IP_1", 0.8656449793155323),
            ("GP_1_1", 0.8656449793155323),
            ("theta_2", 6.283185307179587),
            ("IP_2", 1.3437758555848063),
            ("GP_2_2", 1.3437758555848063),
            ("theta_3", 6.2476281627409715),
            ("IP_3", 0.9405791650996617),
            ("GP_3_3", 0.9405791650996617),
            ("theta_4", 6.127903122194022),
            ("IP_4", 0.0),
            ("theta_5", 6.096864394466671),
            ("IP_5", -0.9),
            ("theta_6", 6.192510223666131),
            ("IP_6", 0.0),
            ("theta_7", 6.154412194314118),
            ("IP_7", -1.0),
            ("theta_8", 6.1991993162055365),
            ("IP_8", 0.0),
            ("theta_9", 6.083000384352472),
            ("IP_9", -1.25),
            ("P_1_1_4", 0.8656449793155323),
            ("twinP_1_1_4", 0.0),
            ("P_2_4_5", 0.3373774752972991),
            ("twinP_2_4_5", 0.0),
            ("P_3_5_6", -0.5626225247027009),
            ("twinP_3_5_6", 0.0),
            ("P_4_3_6", 0.9405791650996617),
            ("twinP_4_3_6", 0.0),
            ("P_5_6_7", 0.37795664039696075),
            ("twinP_5_6_7", 0.0),
            ("P_6_7_8", -0.6220433596030392),
            ("twinP_6_7_8", 0.0),
            ("P_7_8_2", -1.3437758555848063),
            ("twinP_7_8_2", 0.0),
            ("P_8_8_9", 0.7217324959817668),
            ("twinP_8_8_9", 0.0),
            ("P_9_9_4", -0.5282675040182332),
            ("twinP_9_9_4", 0.0),
            ("z_1_1_4", 1.0),
            ("z_2_4_5", 1.0),
            ("z_3_5_6", 1.0),
            ("z_4_3_6", 1.0),
            ("z_5_6_7", 1.0),
            ("z_6_7_8", 1.0),
            ("z_7_8_2", 1.0),
            ("z_8_8_9", 1.0),
            ("z_9_9_4", 1.0),
            ("lincost", 688.1335088379091),
            ("constant", 1.0),
        ]
    )
    objval = 5216.026607747274

    return solution, objval
