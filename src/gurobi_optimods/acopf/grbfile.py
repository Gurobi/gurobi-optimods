import sys
import os
import numpy as np
import time
from myutils import break_exit
import math

def read_configfile(alldata, filename):
    """Function to read configurations for OPF solve from config file"""

    log = alldata['log']

    try:
        log.joint("  Reading configuration file %s\n"%filename)
        f     = open(filename, "r")
        lines = f.readlines()
        f.close()
    except:
        log.raise_exception("Error: Cannot open file %s\n"%filename)

    # Currently Hard-coded
    alldata['usequadcostvar'] = False

    # Default values for options
    casefilename = 'NONE'
    voltsfilename = gvfilename = 'NONE'
    lpfilename = 'grbopf.lp'
    strictcheckvoltagesolution = usevoltsolution = fixcs = skipjabr = cutplane = dodc = False
    useconvexformulation = False
    usemaxdispersion     = False
    use_ef               = False
    substitute_nonconv   = False
    dographics           = False
    dopolar              = False
    doac                 = True
    branchswitching_mip  = False
    branchswitching_comp = False    
    maxdispersion_deg    = 0
    fixtolerance         = 0
    linenum              = 0

    # Read lines of configuration file and save options
    while linenum < len(lines):
        thisline = lines[linenum].split()

        if len(thisline) <= 0:
            linenum += 1
            continue

        if thisline[0] == 'casefilename':
            casefilename = thisline[1]

        elif thisline[0] == 'lpfilename':
            lpfilename = thisline[1]

        elif thisline[0] == 'use_ef':
            use_ef = True

        elif thisline[0] == 'substitute_nonconv':
            substitute_nonconv = True

        elif thisline[0] == 'dopolar':
            dopolar = True

        elif thisline[0] == 'branchswitching_mip':
            branchswitching_mip = True

        elif thisline[0] == 'branchswitching_comp':
            branchswitching_comp = True

        elif thisline[0] == 'usemaxdispersion':
            usemaxdispersion  = True
            maxdispersion_deg = float(thisline[1])

        elif thisline[0] == 'strictcheckvoltagesolution':
            strictcheckvoltagesolution = True
            voltsfilename        = thisline[1]

        elif thisline[0] == 'voltsfilename':
            if len(thisline) < 3:
                log.raise_exception('Error: voltsfilename option requires filename and tolerance\n')
            voltsfilename = thisline[1]
            fixtolerance = float(thisline[2])
            if usevoltsolution:
                log.joint("Cannot use both voltsfilename and voltsolution\n")
                log.raise_exception("Error: Illegal option combination\n")

        elif thisline[0] == 'usevoltsoution':
            if voltsfilename != 'NONE':
                log.joint("Cannot use both voltsfilename and voltsolution\n")
                log.raise_exception("Error: Illegal option combination\n")

            usevoltsolution = True  # forcing solution in casefile


        elif thisline[0] == 'FIXCS':
            fixcs = True

        elif thisline[0] == 'useconvexformulation':
            useconvexformulation = True

        elif thisline[0] == 'gvfileoutput':
            gvfilename = thisline[1]

        elif thisline[0] == 'skipjabr':
            skipjabr = True

        elif thisline[0] == 'cutplane':
            cutplane = True

        elif thisline[0] == 'dodc':
            dodc      = 1
            doac      = 0
            dodcbasic = 0

        elif thisline[0] == 'dodcbasic':
            dodc      = 0
            doac      = 0
            dodcbasic = 1

        elif thisline[0] == 'doac':
            dodc      = 0
            doac      = 1
            dodcbasic = 0

        elif thisline[0] == 'dographics':
            dographics = True

        elif thisline[0] == 'END':
            break

        else:
            log.raise_exception("Error: Illegal input %s\n"%thisline[0])

        linenum += 1

    log.joint("Settings:\n")
    for x in [('casefilename', casefilename), ('lpfilename', lpfilename),
              ('strictcheckvoltagesolution', strictcheckvoltagesolution),
              ('voltsfilename', voltsfilename), ('usevoltsolution', usevoltsolution),
              ('FIXCS', fixcs), ('useconvexformulation', useconvexformulation),
              ('skipjabr', skipjabr), ('cutplane', cutplane), ('dodc',dodc),
              ('doac', doac), ('fixtolerance', fixtolerance), ('use_ef', use_ef),
              ('substitute_nonconv', substitute_nonconv), ('dopolar', dopolar),
              ('usemaxdispersion', usemaxdispersion), ('maxdispersion_deg', maxdispersion_deg),
              ('dographics',dographics), ('branchswitching_mip',branchswitching_mip), ('branchswitching_comp',branchswitching_comp)]:
        alldata[x[0]] = x[1]
        log.joint("  {} {}\n".format(x[0], x[1]))

    if alldata['casefilename'] == 'NONE':
       log.raise_exception('Error: No casefile provided\n')


def grbreadvoltsfile(alldata):
    log = alldata['log']

    numbuses = alldata['numbuses']
    buses = alldata['buses']
    IDtoCountmap = alldata['IDtoCountmap']

    filename = alldata['voltsfilename']
    log.joint("reading volts file " + filename + "\n")

    try:
        f = open(filename, "r")
        lines = f.readlines()
        f.close()
    except:
        log.raise_exception("cannot open file "+filename)


    linenum = 0
    
    while linenum < len(lines):
        thisline = lines[linenum].split()
        if len(thisline) > 0:

            if thisline[0] == 'END':
                break
            elif thisline[0] == 'bus':
                busid = int(thisline[1])
                buscount = IDtoCountmap[busid]

                vm = float(thisline[3])
                vadeg = float(thisline[5])
                varad = vadeg*math.pi/180.
                #print(busid, buscount, vm, vadeg, varad)

                buses[buscount].inputvoltage = 1
                buses[buscount].inputV = vm
                buses[buscount].inputA_rad = varad
            else:
                break_exit('illegal line ' + str(linenum+1))
        linenum += 1
        
    log.joint('Read volts file.\n')

    
       
