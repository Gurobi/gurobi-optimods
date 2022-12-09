import sys
import os
import numpy as np
import time

def read_config(alldata, filename):

    log = alldata['log']

    log.joint("reading config file " + filename + "\n")

    try:
        f = open(filename, "r")
        lines = f.readlines()
        f.close()
    except:
        log.stateandquit("cannot open file", filename)
        sys.exit("failure")

    #hard-coded
    alldata['usequadcostvar'] = False

    casefilename = lpfilename = strictvoltsfilename = 'NONE'
    voltsfilename = gvfilename = 'NONE'
    strictcheckvoltagesolution = usevoltsolution = fixcs = skipjabr = cutplane = dodc = False
    useconvexformulation = False
    usemaxdispersion = False
    maxdispersion_deg = 0
    fixtolerance = 0
    doac = True
    use_ef = False
    substitute_nonconv = False
    dographics = False
    dopolar = False

    linenum = 0
    
    
    while linenum < len(lines):
        thisline = lines[linenum].split()
        if len(thisline) > 0:

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

            elif thisline[0] == 'usemaxdispersion':
                usemaxdispersion = True
                maxdispersion_deg = float(thisline[1])

            elif thisline[0] == 'strictcheckvoltagesolution':
                strictcheckvoltagesolution = True
                strictvoltsfilename = thisline[1]

            elif thisline[0] == 'voltsfilename':
                voltsfilename = thisline[1]
                if usevoltsolution:
                    log.joint("cannot use both voltsfilename and voltsolution\n")
                    log.stateandquit("illegal option\n")

            elif thisline[0] == 'usevoltsoution':
                if voltsfilename != 'NONE':
                    log.joint("cannot use both voltsfilename and voltsolution\n")
                    log.stateandquit("illegal option\n")
                usevoltsolution = True  #forcing solution in casefile

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
                dodc = 1
                doac = 0
                dodcbasic = 0
                
            elif thisline[0] == 'dodcbasic':
                dodcbasic = 1
                doac = 0
                dodc = 0

            elif thisline[0] == 'doac':
                dodcbasic = 0
                doac = 1
                dodc = 0

            elif thisline[0] == 'dographics':
                dographics = True

            elif thisline[0] == 'END':
                break
                
            else:
                sys.exit("illegal input " + thisline[0] + "bye")


        linenum += 1


    for x in [('casefilename', casefilename), ('lpfilename', lpfilename), ('strictcheckvoltagesolution', strictcheckvoltagesolution), ('voltsfilename', voltsfilename), ('usevoltsolution', usevoltsolution), ('FIXCS', fixcs), ('useconvexformulation', useconvexformulation), ('skipjabr', skipjabr), ('cutplane', cutplane), ('dodc',dodc), ('doac', doac), ('fixtolerance', fixtolerance), ('use_ef', use_ef), ('substitute_nonconv', substitute_nonconv), ('dopolar', dopolar), ('usemaxdispersion', usemaxdispersion), ('maxdispersion_deg', maxdispersion_deg), ('dographics',dographics)]:
        alldata[x[0]] = x[1]
        log.joint('{} {}\n'.format(x[0], x[1]))

    if alldata['casefilename'] == 'NONE':
       log.stateandquit('no casefile provided\n')

