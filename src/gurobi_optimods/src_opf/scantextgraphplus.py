import logging
import numpy as np


def scantextgraph(alldata, textfilename):
    logger = logging.getLogger("OpfLogger")
    logger.info("Scanning graph in text file %s\n." % (textfilename))

    f = open(textfilename, "r")
    lines = f.readlines()
    f.close()

    logger.info("Using given lat, lon coordinates.\n")
    N = alldata["numbuses"]
    buses = alldata["buses"]
    nodex = np.zeros(N)
    nodey = np.zeros(N)
    for node in range(N):
        nodex[node] = buses[node + 1].lon
        nodey[node] = buses[node + 1].lat

    N = int(lines[0].split()[1])
    M = int(lines[0].split()[3])
    buses = alldata["buses"]

    trueM = 0
    endbus = {}
    revendbus = {}
    scanned_list_consolidated = {}
    scanned_degrees_consolidated = {}
    scanned_unique_ordered_pairs = {}
    scanned_num_unique = 0
    loud = False
    for i in range(1, M + 1):
        foo = lines[i].split()
        busfrom = int(foo[0])
        busto = int(foo[1])
        small = min(busfrom, busto)
        large = max(busfrom, busto)
        if (small, large) not in scanned_list_consolidated.keys():
            scanned_degrees_consolidated[(small, large)] = 1
            scanned_list_consolidated[(small, large)] = []
            scanned_list_consolidated[(small, large)].append(trueM)
            scanned_unique_ordered_pairs[scanned_num_unique] = (small, large)
            scanned_num_unique += 1
            if loud:
                logger.info(
                    " --> Line %d creates scanned consolidated list for (%d,%d) --> unique ct %d\n"
                    % (trueM, small, large, scanned_num_unique)
                )
        else:
            scanned_degrees_consolidated[(small, large)] += 1
            scanned_list_consolidated[(small, large)].append(trueM)
            if loud:
                logger.info(
                    " --> Appended line %d to scanned consolidated list for (%d,%d)\n"
                    % (trueM, small, large)
                )

        endbus[trueM] = (busfrom, busto)
        revendbus[(busfrom, busto)] = trueM + 1
        revendbus[(busto, busfrom)] = trueM + 1
        trueM += 1
    logger.info("After scanning, number of edges is %d.\n" % trueM)
    return (
        trueM,
        nodex,
        nodey,
        scanned_list_consolidated,
        scanned_degrees_consolidated,
        scanned_unique_ordered_pairs,
        scanned_num_unique,
    )
