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

def graphplot(graphfilename, gvfilename):
    """Description"""#FIXME add more details

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
    break_exit('lmn')

    #trueN, nodex, nodey = scangv('first.gv')
    trueN, N, nodex, nodey = scangv(gvfilename)

    print(len(nodex), len(nodey), trueN, N)

    break_exit('scanned')

    pos = {}
    for j in range(n):
        pos[j] = ( [nodex[j], nodey[j]] )

    G = nx.Graph()

    print(len(adj),m)

    break_exit('lmn2')

    for j in range(truelinect):
        G.add_edge(adj[j][0],adj[j][1])
    print('creating visualization object\n')
    vis = GraphVisualization(G, pos, node_size=1, node_border_width=1, edge_width=1.5)
    print('rendering figure\n')

    xgap     = np.max(nodex) - np.min(nodex)
    ygap     = np.max(nodey) - np.min(nodey)
    myheight = 800
    mywidth  = int(math.ceil(xgap/ygap*myheight))
    print('xgap',xgap,'ygap',ygap,'height',myheight,'width',mywidth)

    fig = vis.create_figure(height=myheight, width=mywidth, showlabel=False)
    #fig.write_image('one.png')
    print('showing figure\n')
    fig.show()
