# https://networkx.org/documentation/stable/_downloads/networkx_reference.pdf

import plotly.graph_objects as go
import plotutils as pu
import psutil
import numpy as np
#import networkx as nx
from grbgraph import *
import math
from graphvisualization import *
from scangvplus import *
from myutils import break_exit

def graphplot(alldata, graphfilename, gvfilename, node_text, mynode_size, mynode_color, myedge_width, myedge_color, myedge_ends, myedge_list_consolidated, myedge_degrees_consolidated, numbranches):
    """Description"""
    #
    # Reads a network in the format created by the graphviz library
    # Then uses the GraphVisualization library to create a plotly figure
    # The figure is then rendered in a browser window
    #
    print('graphfilename', graphfilename, 'gvfilename', gvfilename)
    #graphfilename = 'grbgraph.txt'
    try:
        f     = open(graphfilename, "r")
        lines = f.readlines()
        f.close()
    except:
        sys.exit("failure")

    log = alldata['log']

    n = int(lines[0].split()[1])
    m = int(lines[0].split()[3])

    adj        = {}
    truelinect = 0
    for line in range(1, m+1):
        if lines[line].split()[0] == 'END':
            print('found end of graph file after {} lines\n'.format(truelinect))
            break
        #print(line,int(lines[line].split()[0]), int(lines[line].split()[1]))
        adj[truelinect] = [int(lines[line].split()[0])-1, int(lines[line].split()[1])-1]
        truelinect += 1
    #print('lines',len(lines),'m',m, 'n',n)

    #should check that the next line is 'END'
    #break_exit('lmn')

    trueN, N, trueM, nodex, nodey, thisM, endbus, revendbus, scanned_list_consolidated, scanned_degrees_consolidated, scanned_unique_ordered_pairs, scanned_num_unique = scangv(gvfilename, log)
    #print(len(nodex), len(nodey), trueN, N)

    #endbus is a list of end buses for network branches as ordered in gvfilename
    #revendbus is the reverse index
    #scanned_list_consolidated, degrees_consolidated are dictionaries whose keys are the node pairs (each pair in sorted order) corresponding to
    #lines '--' in the .gv file.  list_consolidated is the array of lines whose ends are that node pair, listing the order in which
    #the lines are found in the .gv file.  degrees_consolidated is the number of such lines (for each sorted node pair found in the file

    if trueM != numbranches:  
        log.joint('Error. Number of branches in .gv file = %d, whereas number of branches in problem = %d\n'%(trueM, numbranches))
        break_exit('error')  #<--- Jarek, I suppose we have to throw an exception here.  Should this happen it indicates a bug.
    if trueM != truelinect:  
        log.joint('Error. Number of branches in .gv file = %d, whereas number of branches line file = %d\n'%(trueM, truelinect))
        break_exit('error')  #<--- Jarek, I suppose we have to throw an exception here.  Should this happen it indicates a bug.


    loud = True
    local_reordered_width = {}
    local_reordered_color = {}
    
    for j in range(scanned_num_unique):
        scannedordpair = scanned_unique_ordered_pairs[j]
        degscanned = scanned_degrees_consolidated[scannedordpair]
        if scannedordpair not in myedge_degrees_consolidated.keys():
            log.joint('Error. Ordered pair %d -> (%d,%d) not in consolidated list\n'%(j, scannedordpair[0],scannedordpair[1]))
            break_exit('error')
        degmine = myedge_degrees_consolidated[scannedordpair]
        if degscanned != degmine:
            log.joint('Error. Ordered pair %d -> (%d,%d) has different degrees %d, %d in scanned vs consolidated lists\n'%(j, scannedordpair[0],scannedordpair[1]), degscanned, degmine)
            break_exit('error')
            
        if False:
            log.joint('scanned ordered pair %d -> (%d,%d) of deg %d\n'%(j, scannedordpair[0],scannedordpair[1], degscanned))

        local_reordered_color[scannedordpair] = {}
        local_reordered_width[scannedordpair] = {}
            
        for k in range(degmine):
            j = myedge_list_consolidated[scannedordpair][k]
            jprime = scanned_list_consolidated[scannedordpair][k]
            color = myedge_color[j]
            width = myedge_width[j]
            
            local_reordered_color[scannedordpair][k] = color  
            local_reordered_width[scannedordpair][k] = width
            
            if loud and width > 1:
                log.joint(' Width %d edge (%d,%d) orderings %d -> %d color %s.\n'%(width, scannedordpair[0], scannedordpair[1], j, jprime, color))


    #break_exit('scanned2')

    pos = {}

    for j in range(n):
        pos[j] = ( [nodex[j], nodey[j]] )

    #G = nx.Graph() #doesn't work too well
    gG = grbGraph(log)

    for j in range(n):
        pos[j] = ( [nodex[j], nodey[j]] )
        gG.addnode(j)
    
    #print(len(adj),m)

    newdeg = {}
    for j in range(scanned_num_unique):
        scannedordpair = scanned_unique_ordered_pairs[j]
        newdeg[scannedordpair[0],scannedordpair[1]] = 0


    #print(newdeg)
    #print(scanned_list_consolidated.keys())

    reordered_width = {}
    reordered_color = {}
    reordered_position = {}

    loud = False
    for j in range(truelinect):
        small = min(adj[j][0],adj[j][1])
        large = max(adj[j][0],adj[j][1])
        #G.add_edge(small, large)
        addcode = gG.addedge(small, large)
        if addcode:
            log.joint('Could not add edge (%d, %d) to graph object.\n'%(small, large))
            break_exit('error')
        fbus = small + 1
        tbus = large + 1
        #print(j+1, 'f,t', fbus, tbus)
        pair = (fbus, tbus)
        deg = newdeg[pair]
        if ((fbus,tbus)) in scanned_list_consolidated.keys():
            if loud and local_reordered_width[pair][deg] > 0:
                log.joint(' pair ind0 %d (%d, %d) local color %s width %d\n'%(j, fbus, tbus, local_reordered_color[pair][deg], local_reordered_width[pair][deg]))
            newdeg[fbus,tbus] += 1
        reordered_width[j] = local_reordered_width[pair][deg]
        reordered_color[j] = local_reordered_color[pair][deg]
        reordered_position[(fbus,tbus,deg)] = j

    gG.getmetrics()
    #break_exit('poo')
    print('Creating visualization object.\n')
    #print(node_text)

    vis = GraphVisualization(gG, pos, node_text, node_size=mynode_size, node_color = mynode_color, node_border_width=1, edge_width = reordered_width, edge_color = reordered_color, edge_map = reordered_position) #myedge_ends) 

    vis.addlog(log)


    #for edge in gG.edges.values():
    #    print(edge[0], edge[1])
    #break_exit('showed')
    
    print('Rendering figure.\n')

    xgap     = np.max(nodex) - np.min(nodex)
    ygap     = np.max(nodey) - np.min(nodey)
    myheight = 800
    mywidth  = int(math.ceil(xgap/ygap*myheight))
    #print('xgap',xgap,'ygap',ygap,'height',myheight,'width',mywidth)

    fig = vis.create_figure(height=myheight, width=mywidth, showlabel=False)
    #fig.write_image('one.png')
    print('Showing figure.\n')
    fig.show()
