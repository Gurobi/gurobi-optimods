import sys
import logging

import gurobipy as gp

from .src_opf.grbcasereader import read_case, build_data_struct

from .src_opf.grbfile import (
    initialize_data_dict,
    read_optimization_settings,
    read_graphics_settings,
    grbread_coords,
    grbread_graphattrs,
)
from .src_opf.grbformulator import construct_and_solve_model
from .src_opf.grbgraphical import plot_solution


def solve_opf_model(settings, case, logfile=""):
    """Construct an OPF model from given data and solve it with Gurobi"""

    if not logfile:
        logfile = "gurobiOPF.log"
    # Initialize output and file handler and start logging
    filehandler = logging.FileHandler(filename=logfile)
    stdouthandler = logging.StreamHandler(stream=sys.stdout)
    handlers = [filehandler, stdouthandler]
    logging.basicConfig(level=logging.INFO, format="%(message)s", handlers=handlers)
    logger = logging.getLogger("OpfLogger")

    # Initilize data dictionary
    alldata = initialize_data_dict(logfile)

    # Read settings file/dict and possibly set case name
    read_optimization_settings(alldata, settings, case)

    # Read case file/dict and populate the alldata dictionary
    read_case(alldata, case)

    solution, objval = construct_and_solve_model(alldata)

    # Close logging handlers
    for handler in handlers:
        handler.close()

    return solution, objval


def plot_opf_solution(settings, case, solution, objval):
    """
    Read the given case and plot a given solution.
    In best case the solution has been computed by the solve_opf_model function
    """

    # Initilize data dictionary
    alldata = initialize_data_dict()

    # Read settings file/dict and possibly set case name
    read_graphics_settings(alldata, settings, case)

    # Read case file/dict and populate the alldata dictionary
    read_case(alldata, case)

    alldata["graphical"] = {}
    alldata["graphical"]["numfeatures"] = 0
    if alldata["graphattrsfilename"] != None:
        grbread_graphattrs(alldata, alldata["graphattrsfilename"])
    if alldata["coordsfilename"] != None:
        grbread_coords(alldata)

    plot_solution(alldata, solution, objval)
