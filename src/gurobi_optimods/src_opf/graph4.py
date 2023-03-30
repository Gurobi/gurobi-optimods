# https://networkx.org/documentation/stable/_downloads/networkx_reference.pdf

import plotly.graph_objects as go
import plotutils as pu
import psutil
import logging
import numpy as np

# import networkx as nx
import math

from .plotlyhandler import *
from .grbgraph import *
from .scangvplus import *
from .myutils import break_exit


def graphplot(
    alldata,
    graphfilename,
    gvfilename,
    myvertex_text,
    myvertex_size,
    myvertex_color,
    myvertex_border_width,
    myedge_width,
    myedge_color,
    myedge_ends,
    myedge_list_consolidated,
    myedge_degrees_consolidated,
    numbranches,
    textlist,
):
    """Description"""
    #
    # When alldata['coordsfilename'] == None, it reads a network in the format created by the graphviz library
    # Then uses the plotlyhandler library to create a plotly figure
    # The figure is then rendered in a browser window
    #
    # print('graphfilename', graphfilename, 'gvfilename', gvfilename)
    # graphfilename = 'grbgraph.txt'
    f = open(graphfilename, "r")
    lines = f.readlines()
    f.close()

    n = int(lines[0].split()[1])
    m = int(lines[0].split()[3])

    adj = {}
    truelinect = 0
    for line in range(1, m + 1):
        if lines[line].split()[0] == "END":
            print("found end of graph file after {} lines\n".format(truelinect))
            break
        # print(line,int(lines[line].split()[0]), int(lines[line].split()[1]))
        adj[truelinect] = [
            int(lines[line].split()[0]) - 1,
            int(lines[line].split()[1]) - 1,
        ]
        truelinect += 1

    # trueN, N, trueM, vertexx, vertexy, thisM, endbus, revendbus, scanned_list_consolidated, scanned_degrees_consolidated, scanned_unique_ordered_pairs, scanned_num_unique = scangv(gvfilename, log)
    # print(len(vertexx), len(vertexy), trueN, N)

    (
        trueM,
        vertexx,
        vertexy,
        scanned_list_consolidated,
        scanned_degrees_consolidated,
        scanned_unique_ordered_pairs,
        scanned_num_unique,
    ) = scangv(alldata, gvfilename, alldata["coordsfilename"] == None)

    # break_exit('scanned')

    # endbus is a list of end buses for network branches as ordered in gvfilename
    # revendbus is the reverse index
    # scanned_list_consolidated, degrees_consolidated are dictionaries whose keys are the vertex pairs (each pair in sorted order) corresponding to
    # lines '--' in the .gv file.  list_consolidated is the array of lines whose ends are that vertex pair, listing the order in which
    # the lines are found in the .gv file.  degrees_consolidated is the number of such lines (for each sorted vertex pair found in the file

    if trueM != numbranches:
        logging.error(
            "Error. Number of branches in .gv file = %d, whereas number of branches in problem = %d."
            % (trueM, numbranches)
        )
        raise ValueError(
            "Number of branches in .gv file = %d, whereas number of branches in problem = %d."
            % (trueM, numbranches)
        )
    if trueM != truelinect:
        logging.error(
            "Error. Number of branches in .gv file = %d, whereas number of branches line file = %d."
            % (trueM, truelinect)
        )
        raise ValueError(
            "Number of branches in .gv file = %d, whereas number of branches line file = %d."
            % (trueM, truelinect)
        )

    loud = False
    local_reordered_width = {}
    local_reordered_color = {}

    for j in range(scanned_num_unique):
        scannedordpair = scanned_unique_ordered_pairs[j]
        degscanned = scanned_degrees_consolidated[scannedordpair]
        if scannedordpair not in myedge_degrees_consolidated.keys():
            logging.error(
                "Error. Ordered pair %d -> (%d,%d) not in consolidated list."
                % (j, scannedordpair[0], scannedordpair[1])
            )
            raise ValueError(
                "Ordered pair %d -> (%d,%d) not in consolidated list."
                % (j, scannedordpair[0], scannedordpair[1])
            )
        degmine = myedge_degrees_consolidated[scannedordpair]
        if degscanned != degmine:
            logging.error(
                "Error. Ordered pair %d -> (%d,%d) has different degrees %d, %d in scanned vs consolidated lists."
                % (j, scannedordpair[0], scannedordpair[1]),
                degscanned,
                degmine,
            )
            raise ValueError(
                "Ordered pair %d -> (%d,%d) has different degrees %d, %d in scanned vs consolidated lists."
                % (j, scannedordpair[0], scannedordpair[1]),
                degscanned,
                degmine,
            )

        if False:  # TODO why is this?
            logging.info(
                "scanned ordered pair %d -> (%d,%d) of deg %d."
                % (j, scannedordpair[0], scannedordpair[1], degscanned)
            )

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
                logging.info(
                    " Width %d edge (%d,%d) orderings %d -> %d color %s."
                    % (width, scannedordpair[0], scannedordpair[1], j, jprime, color)
                )

    # break_exit('scanned2')

    if alldata["coordsfilename"] != None:  # should change this
        vertexx -= np.min(vertexx)
        vertexy -= np.min(vertexy)

        # print('vertexx', vertexx)
        # print('vertexy', vertexy)

        deltax = np.max(vertexx) - np.min(vertexx)  # min should be 0 now
        deltay = np.max(vertexy) - np.min(vertexy)  # min should be 0 now

        ratioxy = deltax / deltay

        scaley = 290
        scalex = 200  # ratio approximately 1.5, mimicking the ratio between one degree of latitute and longitude

        vertexx *= scalex
        vertexy *= scaley

        if False:
            print("vertexx", vertexx)
            print("vertexy", vertexy)

            break_exit("printedvertices")

    pos = {}
    for j in range(n):
        pos[j] = [vertexx[j], vertexy[j]]

    gG = grbGraph()

    for j in range(n):
        gG.addvertex(j)

    # print(len(adj),m)

    newdeg = {}
    for j in range(scanned_num_unique):
        scannedordpair = scanned_unique_ordered_pairs[j]
        newdeg[scannedordpair[0], scannedordpair[1]] = 0

    # print(newdeg)
    # print(scanned_list_consolidated.keys())

    reordered_width = {}
    reordered_color = {}
    reordered_position = {}

    loud = False
    for j in range(truelinect):
        small = min(adj[j][0], adj[j][1])
        large = max(adj[j][0], adj[j][1])
        # G.add_edge(small, large)
        addcode = gG.addedge(small, large)
        if addcode:
            logging.error(
                "Could not add edge (%d, %d) to graph object." % (small, large)
            )
            raise ValueError(
                "Could not add edge (%d, %d) to graph object." % (small, large)
            )
        fbus = small + 1
        tbus = large + 1
        # print(j+1, 'f,t', fbus, tbus)
        pair = (fbus, tbus)
        deg = newdeg[pair]
        if ((fbus, tbus)) in scanned_list_consolidated.keys():
            if loud and local_reordered_width[pair][deg] > 0:
                logging.info(
                    " pair ind0 %d (%d, %d) local color %s width %d."
                    % (
                        j,
                        fbus,
                        tbus,
                        local_reordered_color[pair][deg],
                        local_reordered_width[pair][deg],
                    )
                )
            newdeg[fbus, tbus] += 1
        reordered_width[j] = local_reordered_width[pair][deg]
        reordered_color[j] = local_reordered_color[pair][deg]
        reordered_position[(fbus, tbus, deg)] = j

    gG.getmetrics()
    # break_exit('poo')
    print("Creating visualization object.\n")
    # print(vertex_text)

    PH = plotlyhandler(
        gG,
        pos,
        annotation_list=textlist,
        vertex_size=myvertex_size,
        vertex_color=myvertex_color,
        edge_width=reordered_width,
        edge_color=reordered_color,
        edge_map=reordered_position,
        vertex_text=myvertex_text,
        vertex_border_width=myvertex_border_width,
    )

    # for edge in gG.edges.values():
    #    print(edge[0], edge[1])
    # break_exit('showed')

    logging.info("Rendering figure.")

    xgap = np.max(vertexx) - np.min(vertexx)
    ygap = np.max(vertexy) - np.min(vertexy)

    myheight = 800
    mywidth = int(math.ceil(xgap / ygap * myheight))

    print("xgap", xgap, "ygap", ygap, "height", myheight, "width", mywidth)

    fig = PH.create_figure(height=myheight, width=mywidth, showlabel=False)
    # fig.write_image('one.png')
    logging.info("Showing figure.")
    fig.show()
    # break_exit('showed')
