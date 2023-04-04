import sys
import os
import logging
import re
import numpy as np

from .utils import break_exit


def scangv(alldata, filename, readcoords):
    print("Scanning gv file", filename)

    try:
        f = open(filename, "r")
        lines = f.readlines()
        f.close()
    except:
        return 0, None, None

    N = len(lines)
    # thestrings = np.array(lines[linenum].split() for linenum in range(N))

    criterion = np.array([len(lines[linenum].split()) for linenum in range(N)])
    selection = np.where(criterion == 2)
    buses = alldata["buses"]

    if readcoords:
        N = len(selection[0])
        nodex = np.zeros(N)
        nodey = np.zeros(N)

        for k in range(N):
            i = selection[0][k]
            foo = lines[i].split()
            # print(k, 'i',i, 'split', foo[0], foo[1])
            if len(foo[1]) > 3 and foo[1][:4] == "[pos":
                # print(k, 'i',i, 'split', foo[0], foo[1])
                # print(foo)
                comma = foo[1].rfind(",")
                rightbr = foo[1].rfind("]")
                # print(comma, rightbr)
                node = int(foo[0]) - 1  # base 0
                nodex[node] = float(foo[1][6:comma])
                nodey[node] = float(foo[1][comma + 1 : rightbr - 1])
                # print(node, float(foo[1][6:comma]), float(foo[1][comma+1:rightbr-1]))
    else:
        N = alldata["numbuses"]
        nodex = np.zeros(N)
        nodey = np.zeros(N)
        for node in range(N):
            nodex[node] = buses[node + 1].lon
            nodey[node] = buses[node + 1].lat

    trueM = 0
    N = len(lines)
    endbus = {}
    revendbus = {}
    scanned_list_consolidated = {}
    scanned_degrees_consolidated = {}
    scanned_unique_ordered_pairs = {}
    scanned_num_unique = 0
    loud = False
    for i in range(N):
        foo = lines[i].split()
        # print('k',k,'i',i,len(foo),foo)
        if len(foo) > 2 and foo[1][0] == "-":
            # print(i,len(foo),foo)
            busfrom = int(foo[0])
            busto = int(foo[2])

            small = min(busfrom, busto)
            large = max(busfrom, busto)
            if (small, large) not in scanned_list_consolidated.keys():
                scanned_degrees_consolidated[(small, large)] = 1
                scanned_list_consolidated[(small, large)] = []
                scanned_list_consolidated[(small, large)].append(trueM)
                scanned_unique_ordered_pairs[scanned_num_unique] = (small, large)
                scanned_num_unique += 1
                if loud:
                    logging.info(
                        " --> line %d creates scanned consolidated list for (%d,%d) --> unique ct %d"
                        % (trueM, small, large, scanned_num_unique)
                    )
            else:
                scanned_degrees_consolidated[(small, large)] += 1
                scanned_list_consolidated[(small, large)].append(trueM)
                if loud:
                    logging.info(
                        " --> appended line %d to scanned consolidated list for (%d,%d)"
                        % (trueM, small, large)
                    )

            # print('>>>>>>',busfrom, busto)
            endbus[trueM] = (busfrom, busto)
            revendbus[(busfrom, busto)] = trueM + 1
            revendbus[(busto, busfrom)] = trueM + 1
            trueM += 1
    logging.info("After scanning, number of edges is %d." % trueM)
    # break_exit('scanned')
    return (
        trueM,
        nodex,
        nodey,
        scanned_list_consolidated,
        scanned_degrees_consolidated,
        scanned_unique_ordered_pairs,
        scanned_num_unique,
    )
