import sys
import os
import numpy as np
import time
from log import danoLogger

from grbcasereader import readcase
from myutils import breakexit
from versioner import *
from grbfile import *
from grbgraphical import *
from grbformulator import *

from myconstants import setconstants

if __name__ == '__main__':
    if len(sys.argv) > 3 or len(sys.argv) < 2:
        print ('Usage: grbopf.py file.config [logfile]\n')
        exit(0)
    t0 = time.time()

    mylogfile = "grbopf.log"

    if len(sys.argv) == 3:
        mylogfile = sys.argv[2]

    log = danoLogger(mylogfile)
    stateversion(log)

    alldata = {}
    alldata['log'] = log
    alldata['LP'] = {}

    read_configfile(alldata, sys.argv[1])

    readcase(alldata)

    if alldata['dographics']:
        grbgraphical(alldata)

    lpformulator_ac(alldata)
    breakexit('formulated and solved')

    log.closelog()

def solve_acopf_model(configfile, logfile=""):
    """Construct an ACOPF model from given data and solve it with Gurobi"""

    if not logfile:
        logfile = "gurobiACOPF.log"

    log = danologger(logfile)
    stateversion(log)

    if not isinstance(configfile, str):
        log.stateandquit("Error: Configuration file argument not of type String\n")

    alldata        = {}
    alldata['LP']  = {}
    alldata['log'] = log

    read_configfile(alldata, configfile)
    readcase(alldata)

    #if alldata['dographics']:
    #    grbgraphical(alldata)

    lpformulator_ac(alldata)

    log.closelog()
