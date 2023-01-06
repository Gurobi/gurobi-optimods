# https://networkx.org/documentation/stable/_downloads/networkx_reference.pdf

import plotly.graph_objects as go
import plotutils as pu
import psutil
import numpy as np
import networkx as nx
import math
from graphvisualization import *
from scangvplus import *
from myutils import break_exit

def graphplot(alldata, graphfilename, gvfilename, node_text, mynode_size, mynode_color, myedge_width, myedge_color, myedge_ends):
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
    print('lines',len(lines),'m',m, 'n',n)

    #should check that the next line is 'END'
    #break_exit('lmn')

    trueN, N, nodex, nodey, thisM, endbus = scangv(gvfilename)
    #print(len(nodex), len(nodey), trueN, N)

    #endbus is a list of end buses for network branches as ordered in gvfilename



    #break_exit('scanned')

    pos = {}

    for j in range(n):
        pos[j] = ( [nodex[j], nodey[j]] )

    G = nx.Graph()

    #print(len(adj),m)

    reordered_width = {}
    for j in range(truelinect):
        G.add_edge(adj[j][0],adj[j][1])
        fbus = adj[j][0]+1
        tbus = adj[j][1]+1
        position = myedge_ends[(fbus,tbus)]
        reordered_width[j+1] = 1  #temporary: for debugging only
        #print(j+1, 'f,t', fbus, tbus, 'width', reordered_width[j+1])

    reordered_width[2] = 4
    reordered_width[4] = 4
    #print(reordered_width)

    '''
    for j in range(truelinect):
        position = myedge_ends[endbus[j]]
        print(j+1, 'end',  endbus[j], position)
    '''
        
    break_exit('texted')
    
    print('Creating visualization object.\n')
    #print(node_text)

    vis = GraphVisualization(G, pos, node_text, node_size=mynode_size, node_color = mynode_color, node_border_width=1, edge_width = myedge_width, edge_color = myedge_color, edge_map = myedge_ends) 

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
