import sys

from log import danoLogger
from grbcasereader import readcase
from myutils import breakexit
from grbfile import read_configfile
from grbgraphical import grbgraphical
from grbformulator import lpformulator_ac

if __name__ == '__main__':
    if len(sys.argv) > 3 or len(sys.argv) < 2:
        print ('Usage: grbopf.py file.config [logfile]\n')
        exit(0)

    mylogfile = "grbopf.log"

    if len(sys.argv) == 3:
        mylogfile = sys.argv[2]

    try:
        log = danoLogger(mylogfile)

        alldata        = {}
        alldata['LP']  = {}
        alldata['log'] = log

        read_configfile(alldata, sys.argv[1])

        readcase(alldata)

        if alldata['dographics']:
            grbgraphical(alldata)

        lpformulator_ac(alldata)
        breakexit("formulated and solved")
        log.closelog()

    except gp.GurobiError as e:
        print("Error in Gurobi: Error code %s: %s"%(e.errno, e))

    except Exception as e:
        print(e)


def solve_acopf_model(configfile, logfile = ""):
    """Construct an ACOPF model from given data and solve it with Gurobi"""
    try:
        if not logfile:
            logfile = "gurobiACOPF.log"

        # Create log object
        log = danologger(logfile)

        if not isinstance(configfile, str):
            log.stateandquit("Error: Configuration file argument not of type String\n")

        # Initialize data
        alldata        = {}
        alldata['LP']  = {}
        alldata['log'] = log

        # Read configuration file
        read_configfile(alldata, configfile)

        # Read case file holding OPF network data
        readcase(alldata)

        if alldata['dographics']:
            grbgraphical(alldata)

        # Construct model from collected data and optimize it
        lpformulator_ac(alldata)

        log.closelog()

    except gp.GurobiError as e:
        print("Error in Gurobi: Error code %s: %s"%(e.errno, e))

    except:
        print("")  
