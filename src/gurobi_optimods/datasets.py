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


def load_simpleopfsettings():
    conf = str(DATA_FILE_DIR) + "/opf/kit_simpleopf.conf"
    return conf


def load_acopfsettings():
    conf = str(DATA_FILE_DIR) + "/opf/kit_acopf.conf"
    return conf


def load_dcopfsettings():
    conf = str(DATA_FILE_DIR) + "/opf/kit_dcopf.conf"
    return conf


def load_ivopfsettings():
    conf = str(DATA_FILE_DIR) + "/opf/kit_ivopf.conf"
    return conf


def load_case9opf():
    case = str(DATA_FILE_DIR) + "/opf/case9.m"
    return case


def load_opfdictsettings():
    conf = {"branchswitching_mip": True, "doac": True}
    return conf


def load_opfdictgraphicssettings():
    conf = {"branchswitching_mip": True}
    return conf


def load_opfgraphicssettings():
    conf = str(DATA_FILE_DIR) + "/opf/kit_graphics.conf"
    return conf


def load_opfdictcase():
    casefile_dict = {
        "refbus": 1,
        "baseMVA": 100.0,
        "buses": {
            1: {
                "count": 1,
                "nodeID": 1,
                "nodetype": 3,
                "Pd": 0.0,
                "Qd": 0.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "Vbase": 345.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
                "lnum": 18,
            },
            2: {
                "count": 2,
                "nodeID": 2,
                "nodetype": 2,
                "Pd": 0.0,
                "Qd": 0.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "Vbase": 345.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
                "lnum": 19,
            },
            3: {
                "count": 3,
                "nodeID": 3,
                "nodetype": 2,
                "Pd": 0.0,
                "Qd": 0.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "Vbase": 345.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
                "lnum": 20,
            },
            4: {
                "count": 4,
                "nodeID": 4,
                "nodetype": 1,
                "Pd": 0.0,
                "Qd": 0.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "Vbase": 345.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
                "lnum": 21,
            },
            5: {
                "count": 5,
                "nodeID": 5,
                "nodetype": 1,
                "Pd": 90.0,
                "Qd": 30.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "Vbase": 345.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
                "lnum": 22,
            },
            6: {
                "count": 6,
                "nodeID": 6,
                "nodetype": 1,
                "Pd": 0.0,
                "Qd": 0.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "Vbase": 345.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
                "lnum": 23,
            },
            7: {
                "count": 7,
                "nodeID": 7,
                "nodetype": 1,
                "Pd": 100.0,
                "Qd": 35.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "Vbase": 345.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
                "lnum": 24,
            },
            8: {
                "count": 8,
                "nodeID": 8,
                "nodetype": 1,
                "Pd": 0.0,
                "Qd": 0.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "Vbase": 345.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
                "lnum": 25,
            },
            9: {
                "count": 9,
                "nodeID": 9,
                "nodetype": 1,
                "Pd": 125.0,
                "Qd": 50.0,
                "Gs": 0.0,
                "Bs": 0.0,
                "Vbase": 345.0,
                "Vmax": 1.1,
                "Vmin": 0.9,
                "lnum": 26,
            },
        },
        "slackbus": 1,
        "generators": {
            1: {
                "gencount1": 1,
                "nodeID": 1,
                "Pg": 0.0,
                "Qg": 0.0,
                "status": 1,
                "Pmax": 250.0,
                "Pmin": 10.0,
                "Qmax": 300.0,
                "Qmin": -300.0,
                "lnum": 32,
            },
            2: {
                "gencount1": 2,
                "nodeID": 2,
                "Pg": 163.0,
                "Qg": 0.0,
                "status": 1,
                "Pmax": 300.0,
                "Pmin": 10.0,
                "Qmax": 300.0,
                "Qmin": -300.0,
                "lnum": 33,
            },
            3: {
                "gencount1": 3,
                "nodeID": 3,
                "Pg": 85.0,
                "Qg": 0.0,
                "status": 1,
                "Pmax": 270.0,
                "Pmin": 10.0,
                "Qmax": 300.0,
                "Qmin": -300.0,
                "lnum": 34,
            },
        },
        "branches": {
            1: {
                "branchcount1": 1,
                "f": 1,
                "t": 4,
                "r": 0.0,
                "x": 0.0576,
                "bc": 0.0,
                "rateA": 250.0,
                "rateB": 250.0,
                "rateC": 250.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "minangle": -360.0,
                "maxangle": 360.0,
                "lnum": 40,
            },
            2: {
                "branchcount1": 2,
                "f": 4,
                "t": 5,
                "r": 0.017,
                "x": 0.092,
                "bc": 0.158,
                "rateA": 250.0,
                "rateB": 250.0,
                "rateC": 250.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "minangle": -360.0,
                "maxangle": 360.0,
                "lnum": 41,
            },
            3: {
                "branchcount1": 3,
                "f": 5,
                "t": 6,
                "r": 0.039,
                "x": 0.17,
                "bc": 0.358,
                "rateA": 150.0,
                "rateB": 150.0,
                "rateC": 150.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "minangle": -360.0,
                "maxangle": 360.0,
                "lnum": 42,
            },
            4: {
                "branchcount1": 4,
                "f": 3,
                "t": 6,
                "r": 0.0,
                "x": 0.0586,
                "bc": 0.0,
                "rateA": 300.0,
                "rateB": 300.0,
                "rateC": 300.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "minangle": -360.0,
                "maxangle": 360.0,
                "lnum": 43,
            },
            5: {
                "branchcount1": 5,
                "f": 6,
                "t": 7,
                "r": 0.0119,
                "x": 0.1008,
                "bc": 0.209,
                "rateA": 150.0,
                "rateB": 150.0,
                "rateC": 150.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "minangle": -360.0,
                "maxangle": 360.0,
                "lnum": 44,
            },
            6: {
                "branchcount1": 6,
                "f": 7,
                "t": 8,
                "r": 0.0085,
                "x": 0.072,
                "bc": 0.149,
                "rateA": 250.0,
                "rateB": 250.0,
                "rateC": 250.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "minangle": -360.0,
                "maxangle": 360.0,
                "lnum": 45,
            },
            7: {
                "branchcount1": 7,
                "f": 8,
                "t": 2,
                "r": 0.0,
                "x": 0.0625,
                "bc": 0.0,
                "rateA": 250.0,
                "rateB": 250.0,
                "rateC": 250.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "minangle": -360.0,
                "maxangle": 360.0,
                "lnum": 46,
            },
            8: {
                "branchcount1": 8,
                "f": 8,
                "t": 9,
                "r": 0.032,
                "x": 0.161,
                "bc": 0.306,
                "rateA": 250.0,
                "rateB": 250.0,
                "rateC": 250.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "minangle": -360.0,
                "maxangle": 360.0,
                "lnum": 47,
            },
            9: {
                "branchcount1": 9,
                "f": 9,
                "t": 4,
                "r": 0.01,
                "x": 0.085,
                "bc": 0.176,
                "rateA": 250.0,
                "rateB": 250.0,
                "rateC": 250.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "minangle": -360.0,
                "maxangle": 360.0,
                "lnum": 48,
            },
        },
        "generator_cost_structure": {
            1: {"costtype": 2, "degree": 2, "costvector": [1100.0, 500.0, 150.0]},
            2: {
                "costtype": 2,
                "degree": 2,
                "costvector": [850.0000000000001, 120.0, 600.0],
            },
            3: {"costtype": 2, "degree": 2, "costvector": [1225.0, 100.0, 335.0]},
        },
        "generator_cost_count": 3,
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
