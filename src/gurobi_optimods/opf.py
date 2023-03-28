import sys
import logging

import gurobipy as gp

from .src_opf.grbcasereader import read_case

from .src_opf.grbfile import read_configfile, grbread_coords, gbread_graphattrs
from .src_opf.grbgraphical import grbgraphical
from .src_opf.grbformulator_ac import lpformulator_ac
from .src_opf.grbformulator_dc import lpformulator_dc
from .src_opf.grbformulator_iv import lpformulator_iv


def solve_opf_model(configfile, casefile, logfile=""):
    """Construct an ACOPF model from given data and solve it with Gurobi"""

    if not logfile:
        logfile = "gurobiOPF.log"

    filehandler = logging.FileHandler(filename=logfile)
    stdouthandler = logging.StreamHandler(stream=sys.stdout)

    handlers = [filehandler, stdouthandler]

    logging.basicConfig(level=logging.INFO, format="%(message)s", handlers=handlers)

    alldata = {}
    alldata["LP"] = {}
    alldata["MIP"] = {}
    alldata["logfile"] = logfile

    # Read configuration file
    read_configfile(alldata, configfile, casefile)

    # Read case file holding OPF network data
    read_case(alldata)

    if alldata["dographics"]:
        alldata["graphical"] = {}
        alldata["graphical"]["numfeatures"] = 0
        if alldata["graphattrsfilename"] != None:
            gbread_graphattrs(alldata, alldata["graphattrsfilename"])
        if alldata["coordsfilename"] != None:
            grbread_coords(alldata)

    # Construct model from collected data and optimize it
    if alldata["doac"]:
        solution, objval = lpformulator_ac(alldata)
    elif alldata["dodc"]:
        solution, objval = lpformulator_dc(alldata)
    elif alldata["doiv"]:
        lpformulator_iv(alldata)
    elif alldata["doslp_polar"]:
        alldata["doac"] = True
        solution, objval = lpformulator_ac(alldata)

    return solution, objval
