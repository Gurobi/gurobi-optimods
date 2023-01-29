import sys

import gurobipy as gp
from log import Logger
#from opfexception import OPFException  #disable for now
from grbcasereader import read_case
from myutils import break_exit

from grbfile import read_configfile, grbread_coords, gbread_graphattrs
from grbgraphical import grbgraphical
from grbformulator_ac import lpformulator_ac
from grbformulator_dc import lpformulator_dc

if __name__ == '__main__':
    if len(sys.argv) > 3 or len(sys.argv) < 2:
        print ('Usage: grbopf.py file.config [logfile]\n')
        exit(0)

    mylogfile = "grbopf.log"

    if len(sys.argv) == 3:
        mylogfile = sys.argv[2]

    log = Logger(mylogfile)

    alldata        = {}
    alldata['LP']  = {}
    alldata['MIP']  = {}        
    alldata['log'] = log

    read_configfile(alldata, sys.argv[1])

    read_case(alldata)

    if alldata['dographics']:
        alldata['graphical'] = {}
        alldata['graphical']['numfeatures'] = 0
        if alldata['graphattrsfilename'] != None:
            gbread_graphattrs(alldata, alldata['graphattrsfilename'])
        if alldata['coordsfilename'] != None:
            grbread_coords(alldata)

    if alldata['doac']:
            lpformulator_ac(alldata)
    else:
        lpformulator_dc(alldata)
            
    break_exit("formulated and solved")
    log.close_log()

    '''
    #Jarek, I disabled this temporarily because the exception handling was giving me problems
    try:
        log = Logger(mylogfile)

        alldata        = {}
        alldata['LP']  = {}
        alldata['MIP']  = {}        
        alldata['log'] = log

        read_configfile(alldata, sys.argv[1])

        read_case(alldata)

        if alldata['dographics'] and alldata['coordsfilename'] != None:
            grbread_coords(alldata)

        if alldata['doac']:
            lpformulator_ac(alldata)
        else:
            lpformulator_dc(alldata)
            
        break_exit("formulated and solved")
        log.close_log()


    except gp.GurobiError as e:
        print("Error in Gurobi: Error code %s: %s"%(e.errno, e))


    except OPFException as e:
        print(str(e))
    '''
    
def solve_acopf_model(configfile, logfile = ""):
    """Construct an ACOPF model from given data and solve it with Gurobi"""

    try:
        if not logfile:
            logfile = "gurobiACOPF.log"

        # Create log object
        log = Logger(logfile)

        if not isinstance(configfile, str):
            log.raise_exception("Error: Configuration file argument not of type String\n")


        # Read configuration file
        read_configfile(alldata, configfile)

        # Read case file holding OPF network data
        read_case(alldata)

        if alldata['dographics'] and alldata['coordsfilename'] != None:
            grbread_coords(alldata)

        # Construct model from collected data and optimize it
        lpformulator_ac(alldata)

        log.close_log()

    except gp.GurobiError as e:
        print("Error in Gurobi: Error code %s: %s"%(e.errno, e))

    except OPFException as e:
        print(str(e))
