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
from .src_opf.grbformulator import construct_and_solve_model


def solve_opf_model(settings, case, logfile=""):
    """Construct an ACOPF model from given data and solve it with Gurobi"""

    if not logfile:
        logfile = "gurobiOPF.log"
    # Initialize output and file handler and start logging
    filehandler = logging.FileHandler(filename=logfile)
    stdouthandler = logging.StreamHandler(stream=sys.stdout)
    handlers = [filehandler, stdouthandler]
    logging.basicConfig(level=logging.INFO, format="%(message)s", handlers=handlers)

    # Initilize data dictionary
    alldata = initialize_data_dict(logfile)

    # Read settings file/dict and possibly set case name
    read_settings(alldata, settings, case)

    # Read case file/dict and populate the alldata dictionary
    read_case(alldata, case)

    # Special settings for graphics
    if alldata["dographics"]:
        alldata["graphical"] = {}
        alldata["graphical"]["numfeatures"] = 0
        if alldata["graphattrsfilename"] != None:
            grbread_graphattrs(alldata, alldata["graphattrsfilename"])
        if alldata["coordsfilename"] != None:
            grbread_coords(alldata)

    solution, objval = construct_and_solve_model(alldata)

    # Close logging handlers
    for handler in handlers:
        handler.close()

    return solution, objval
