import sys
import logging

import gurobipy as gp

from .src_opf.grbcasereader import read_case, build_data_struct

from .src_opf.grbfile import (
    initialize_data_dict,
    read_settings,
    grbread_coords,
    grbread_graphattrs,
)
from .src_opf.grbgraphical import grbgraphical
from .src_opf.grbformulator_ac import lpformulator_ac
from .src_opf.grbformulator_dc import lpformulator_dc
from .src_opf.grbformulator_iv import lpformulator_iv


def solve_opf_model(settings, case, logfile=""):
    """Construct an ACOPF model from given data and solve it with Gurobi"""

    if not logfile:
        logfile = "gurobiOPF.log"
    # Initialize output and file handler and start logging
    filehandler = logging.FileHandler(filename=logfile)
    stdouthandler = logging.StreamHandler(stream=sys.stdout)
    handlers = [filehandler, stdouthandler]
    logging.basicConfig(level=logging.INFO, format="%(message)s", handlers=handlers)

    alldata = initialize_data_dict(logfile)

    # Read configuration file/dict and possibly set case name
    read_settings(alldata, settings, case)

    # Read case and populate the alldata dictionary
    read_case(alldata, case)

    # Special settings for graphics
    if alldata["dographics"]:
        alldata["graphical"] = {}
        alldata["graphical"]["numfeatures"] = 0
        if alldata["graphattrsfilename"] != None:
            grbread_graphattrs(alldata, alldata["graphattrsfilename"])
        if alldata["coordsfilename"] != None:
            grbread_coords(alldata)

    # Construct model from collected data and optimize it
    if alldata["doac"]:
        solution, objval = lpformulator_ac(alldata)
    elif alldata["dodc"]:
        solution, objval = lpformulator_dc(alldata)
    elif alldata["doiv"]:
        solution, objval = lpformulator_iv(alldata)
    elif alldata["doslp_polar"]:
        alldata["doac"] = True
        solution, objval = lpformulator_ac(alldata)

    # Close logging handlers
    for handler in handlers:
        handler.close()

    return solution, objval
