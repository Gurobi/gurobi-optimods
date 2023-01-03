import sys
import math
import time
from graph4 import *
#from log import Logger
from gurobipy import *
from myutils import break_exit

def grbgraphical(alldata):
    log         = alldata['log']
    buses       = alldata['buses']
    numbuses    = alldata['numbuses']
    numbranches = alldata['numbranches']
    gvfilename  = 'grbgraphical.gv'
    txtfilename = 'newgraph.txt'
    log.joint("Graphical layout, 2\n")

    try:
        f = open(gvfilename, "w")
        log.joint("Writing to gv file %s\n"%gvfilename)
    except:
        log.raise_exception("Error: Cannot open file %s\n"%gvfilename)

    try:
        g = open(txtfilename, "w")
        log.joint("Writing to txt file %s\n"%txtfilename)
    except:
        log.raise_exception("Error: Cannot open file %s\n"%txtfilename)

    f.write("graph {\n")
    f.write('node [color=black, height=0, label=\"\\N\", shape=point, width=0];\n')
    g.write('N '+str(numbuses) + ' M ' + str(numbranches) + '\n')
    for bus in alldata['buses'].values():
        #f.write("     " + str(bus.nodeID)+";\n")
        f.write("     " + str(bus.count)+";\n")
    for branch in alldata['branches'].values():
        f.write("     " + str(branch.count_f)+" -- " + str(branch.count_t)+";\n")
        g.write(' ' + str(branch.count_f)+ ' ' + str(branch.count_t)+ '\n')

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

    #break_exit('graph,1')

    Vmagviol = alldata['violation']['Vmagviol']
    IPviol = alldata['violation']['IPviol']
    IQviol = alldata['violation']['IQviol']
    branchlimitviol = alldata['violation']['branchlimit']

    node_text = {}
    mynode_size = {}
    mynode_color = {}


    for j in range(1,numbuses+1):
        bus = buses[j]
        #node_text[j-1] = 'Bus ' + str(j) + ' Vmagviol: '+ str(Vmagviol[bus]) + ' Pviol: '+ str(IPviol[bus]) + ' Qviol: '+ str(IQviol[bus])

        node_text[j-1] = 'Bus %d Vmagviol: %.3e Pviol %.3e Qviol %.3e'%(j, Vmagviol[bus], IPviol[bus],IQviol[bus])
        mynode_size[j-1] = 1
        mynode_color[j-1] = 'black'
        if abs(Vmagviol[bus]) > 1e-3 or abs(IPviol[bus]) > 1e-2 or abs(IQviol[bus]) > 1e-2:
            mynode_size[j-1] = 15
            mynode_color[j-1] = 'red'

    myedge_width = {}
    myedge_ends = {}
    for j in range(1,numbranches+1):
        branch = alldata['branches'][j]
        myedge_ends[(branch.count_f, branch.count_t)] = j
        myedge_ends[(branch.count_t, branch.count_f)] = j        
        myedge_width[j] = 2
        if abs(branchlimitviol[branch]) > 1e-3:
            myedge_width[j] = 4
            
            
        
    graphplot(alldata, txtfilename, firstgvfile, node_text, mynode_size, mynode_color, myedge_width, myedge_ends)
    break_exit('graph,2')
