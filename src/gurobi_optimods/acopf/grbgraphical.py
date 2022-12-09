import sys
import math
import time
from graph4 import *
#from log import Logger
from gurobipy import *
from myutils import breakexit

def grbgraphical(alldata):
    log         = alldata['log']
    numbuses    = alldata['numbuses']
    numbranches = alldata['numbranches']
    gvfilename  = 'grbgraphical.gv'
    txtfilename = 'newgraph.txt'
    log.joint("Graphical layout, 2\n")

    try:
        f = open(gvfilename, "w")
        log.joint("Writing to gv file " + gvfilename + "\n")
    except:
        log.stateandquit("Error: Cannot open file " + gvfilename)

    try:
        g = open(txtfilename, "w")
        log.joint("Writing to txt file " + txtfilename + "\n")
    except:
        log.stateandquit("Error: Cannot open file " + txtfilename)

    f.write("graph {\n")
    f.write('node [color=black, height=0, label=\"\\N\", shape=point, width=0];\n')
    g.write('N '+str(numbuses) + ' M ' + str(numbranches) + '\n')
    for bus in alldata['buses'].values():
        #f.write("     " + str(bus.nodeID)+";\n")
        f.write("     " + str(bus.count)+";\n")
    for branch in alldata['branches'].values():
        f.write("     " + str(branch.id_f)+" -- " + str(branch.id_t)+";\n")
        g.write(' ' + str(branch.id_f)+ ' ' + str(branch.id_t)+ '\n')

    f.write("}\n")
    f.close()
    g.write('END\n')
    g.close()

    scale       = 10
    firstgvfile = 'first.gv'
    sfdpcommand = 'sfdp -Goverlap_scaling='+str(scale)+' -o '+firstgvfile + ' ' + gvfilename
    log.joint('sfdp command: ' + sfdpcommand + '\n')
    system(sfdpcommand)

    '''
    jpgfile = 'second.jpg'
    neatocommand = 'neato -Tjpeg -n -o '+ jpgfile + ' ' + firstgvfile
    log.joint('neato command: ' + neatocommand + '\n')
    system(neatocommand)
    '''

    breakexit('graph,1')

    graphplot(txtfilename, firstgvfile)

    breakexit('graph,2')
