#!/usr/bin/python
import sys
import os
import numpy as np
import time
from log import danoLogger

import grbcasereader
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

    read_config(alldata,sys.argv[1])

    readcode = grbcasereader.readcase(alldata)
    if readcode:
        sys.exit('cannot read inputfile. \n')

    if alldata['dographics']:
        grbgraphical(alldata)
    
    lpformulator_ac(alldata)
    breakexit('formulated and solved')

    log.closelog()
    
