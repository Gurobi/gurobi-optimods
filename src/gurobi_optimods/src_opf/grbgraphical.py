import logging
import math
import numpy as np

from .grbgraph import Grbgraph
from .plotlyhandler import Plotlyhandler


def generate_solution_figure(alldata, solution, objval):
    """
    Plots a given OPF solution.
    Saves necessary data from solution to alldata dictionary
    and calls the main plotting function.


    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    solution : dictionary
        Dictionary holding an OPF solution
    objval: float
        Objective solution value of a previously solved OPF

    Returns
    -------
    plotly.graph_objects.Figure
        A plotly figure objects which can be displayed via the show() function,
        see https://plotly.com/python-api-reference/generated/plotly.graph_objects.Figure.html
    """

    logger = logging.getLogger("OpfLogger")
    logger.info("Plotting solution with value %.3e. Coordinates given." % objval)
    numbranches = alldata["numbranches"]
    numgens = alldata["numgens"]

    # Save branch and generator values
    zholder = np.zeros(numbranches)
    alldata["MIP"]["zholder"] = zholder
    gholder = np.zeros(alldata["numgens"])
    alldata["MIP"]["gholder"] = gholder
    zholder = alldata["MIP"]["zholder"]
    gholder = alldata["MIP"]["gholder"]

    numzeros = 0
    for j in range(1, 1 + numbranches):
        branch = solution["branch"][j]
        zholder[j - 1] = branch["switching"]
        if branch["switching"] < 0.5:
            numzeros += 1

    for j in range(1, 1 + numgens):
        gen = solution["gen"][j]
        gholder[j - 1] = gen["Pg"]

    textlist = []
    textlist.append("OBJ: %10.2f" % (objval))
    textlist.append("Lines off: %d" % (numzeros))
    return grbgraphical(alldata, "branchswitching", textlist)


def grbgraphical(alldata, plottype, textlist):
    """
    Plots a given OPF solution. Depending on the plottype,
    the edge and node colors are different and we print different
    information

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    plottype : string
        Type of plot to show. Currently, "branchswitching" or "violation"
    textlist: list
        List of additional strings printed on the plot

    Returns
    -------
    plotly.graph_objects.Figure
        A plotly figure objects which can be displayed via the show() function,
        see https://plotly.com/python-api-reference/generated/plotly.graph_objects.Figure.html
    """

    logger = logging.getLogger("OpfLogger")
    buses = alldata["buses"]
    numbuses = alldata["numbuses"]
    branches = alldata["branches"]
    numbranches = alldata["numbranches"]
    gens = alldata["gens"]
    IDtoCountmap = alldata["IDtoCountmap"]
    graph_dict = {}
    logger.info("Graphical layout, coordinates given.\n")

    graph_dict["N"] = numbuses
    graph_dict["M"] = numbranches
    counter = 0
    for branch in alldata["branches"].values():
        graph_dict[counter] = (branch.count_f, branch.count_t)
        counter += 1

    node_text = {}
    mynode_size = {}
    mynode_color = {}
    mynode_border_width = {}
    myedge_width = {}
    myedge_color = {}

    # default actions
    for j in range(1, numbuses + 1):
        bus = buses[j]
        mynode_size[j - 1] = 1
        mynode_color[j - 1] = "black"
        mynode_border_width[j - 1] = 1

    for j in range(1, numbranches + 1):
        branch = branches[j]
        myedge_width[j] = 1
        myedge_color[j] = "black"

    if plottype == "violation":
        Vmagviol = alldata["violation"]["Vmagviol"]
        IPviol = alldata["violation"]["IPviol"]
        IQviol = alldata["violation"]["IQviol"]
        branchlimitviol = alldata["violation"]["branchlimit"]

        for j in range(1, numbuses + 1):
            bus = buses[j]
            # node_text[j-1] = 'Bus ' + str(j) + ' Vmagviol: '+ str(Vmagviol[bus]) + ' Pviol: '+ str(IPviol[bus]) + ' Qviol: '+ str(IQviol[bus])

            node_text[j - 1] = "Bus %d Vmagviol: %.3e Pviol %.3e Qviol %.3e" % (
                j,
                Vmagviol[bus],
                IPviol[bus],
                IQviol[bus],
            )
            if (
                abs(Vmagviol[bus]) > 1e-3
                or abs(IPviol[bus]) > 1e-2
                or abs(IQviol[bus]) > 1e-2
            ):
                mynode_size[j - 1] = 15
                mynode_color[j - 1] = "red"

        for j in range(1, numbranches + 1):
            branch = alldata["branches"][j]
            if abs(branchlimitviol[branch]) > 1e-3:
                myedge_width[j] = 8
                myedge_color[j] = "red"

    elif plottype == "branchswitching":

        loud = False

        gholder = alldata["MIP"]["gholder"]
        for j in range(1, numbuses + 1):
            bus = buses[j]
            sumPgen = 0

            largegen = False

            for k in range(len(bus.genidsbycount)):
                gen = gens[bus.genidsbycount[k]]
                sumPgen += gholder[gen.count - 1]
                if False:
                    logger.info(
                        " bus %d gen %d produces %f\n"
                        % (j, gen.count, gholder[gen.count - 1])
                    )
                # break_exit('gen in grbgraphical')

            """
            if sumPgen > 500:
                largegen = True
            elif sumPgen > 150:
                largegen = True
            elif sumPgen > 100:
                largegen = True
            if False and sumPgen > 0:
                print('bus',j,'sumP', sumPgen, 'color',mynode_color[j-1])
                break_exit('huhhhh')
            """

            largegen, mynode_size[j - 1], mynode_color[j - 1] = grbgetgraphattr(
                alldata, sumPgen
            )

            Pload = alldata["baseMVA"] * bus.Pd
            node_text[j - 1] = "Bus %d   Gen %7.2f  Load %7.2f" % (
                bus.nodeID,
                sumPgen,
                Pload,
            )

            if Pload > 50 and largegen == False:
                mynode_size[j - 1] = 8
                mynode_color[j - 1] = "blue"

        # break_exit('bus examination')

        zholder = alldata["MIP"]["zholder"]
        for j in range(1, 1 + numbranches):
            branch = branches[j]
            f = branch.f
            t = branch.t
            count_of_f = IDtoCountmap[f]
            count_of_t = IDtoCountmap[t]

            if zholder[j - 1] < 0.5:  # turned off
                if loud:
                    logger.info(
                        "branch %d (%d, %d) has small x\n"
                        % (j, branch.count_f, branch.count_t)
                    )
                myedge_width[j] = 10
                myedge_color[j] = "red"

                if mynode_size[count_of_f - 1] < 10:  # a lot of hard-coded values...
                    mynode_size[count_of_f - 1] = 7
                    mynode_color[count_of_f - 1] = "slategrey"

                if mynode_size[count_of_t - 1] < 10:  # a lot of hard-coded values...
                    mynode_size[count_of_t - 1] = 7
                    mynode_color[count_of_t - 1] = "slategrey"

    myedge_ends = {}
    myedge_list_consolidated = {}
    myedge_degrees_consolidated = {}

    loud = False
    for j in range(1, numbranches + 1):
        branch = alldata["branches"][j]
        myedge_ends[(branch.count_f, branch.count_t)] = j
        myedge_ends[(branch.count_t, branch.count_f)] = j
        small = min(branch.count_f, branch.count_t)
        large = max(branch.count_f, branch.count_t)
        if (small, large) not in myedge_list_consolidated.keys():
            myedge_degrees_consolidated[(small, large)] = 1
            myedge_list_consolidated[(small, large)] = []
            myedge_list_consolidated[(small, large)].append(j)
            if loud:
                logger.info(
                    " --> line %d color %s creates my consolidated list for (%d,%d)\n"
                    % (j, myedge_color[j], small, large)
                )
        else:
            myedge_degrees_consolidated[(small, large)] += 1
            myedge_list_consolidated[(small, large)].append(j)
            if loud:
                logger.info(
                    " --> appended line %d color %s to my consolidated list for (%d,%d)\n"
                    % (j, myedge_color[j], small, large)
                )

    # break_exit('cons')

    if False:
        for j in range(1, numbuses + 1):
            if mynode_size[j - 1] > 1:
                print(
                    "pre graphplot, bus",
                    j,
                    "size",
                    mynode_size[j - 1],
                    "color",
                    mynode_color[j - 1],
                )

    return graphplot(
        alldata,
        graph_dict,
        node_text,
        mynode_size,
        mynode_color,
        mynode_border_width,
        myedge_width,
        myedge_color,
        myedge_ends,
        myedge_list_consolidated,
        myedge_degrees_consolidated,
        textlist,
    )


def graphplot(
    alldata,
    graph_dict,
    myvertex_text,
    myvertex_size,
    myvertex_color,
    myvertex_border_width,
    myedge_width,
    myedge_color,
    myedge_ends,
    myedge_list_consolidated,
    myedge_degrees_consolidated,
    textlist,
):
    """
    Plots an OPF network out of user given solution/violation data which has been
    translated into necessary dictionaries in function grbgraphical.py:grbgraphical.
    Uses lat, lon coordinates as in coordinates dictionary,
    together with the plotlyhandler library to create a plotly figure.
    The figure is then rendered in a browser window.

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    graph_dict : dictionary
        Dictionary holding graph data
    myvertex_text : dictionary
        Dictionary holding the texts when hovering over a node
    myvertex_size : dictionary
        Dictionary holding the size of each node. The size depends on generated power
    myvertex_color : dictionary
        Dictionary holding the color of each node. The color depends on generated power
    myvertex_border_width : dictionary
        Dictionary holding the thickness of the border of each node
    myedge_width : dictionary
        Dictionary holding the thickness of each network edge
    myedge_color : dictionary
        Dictionary holding the color of each network edge. The color depends on violation and whether the edge is turned on/off
    myedge_ends : dictionary # TODO-Dan do we need it?
        Dictionary holding each edge in both directions, e.g., (1,2) and (2,1)
    myedge_list_consolidated : dictionary # TODO-Dan What exactly is this?
        Dictionary holding possible multi-edges
    myedge_degrees_consolidated : dictionary
        Dictionary holding the degree of each edge
    textlist : list
        List with additional strings printed on the plot

    Returns
    -------
    plotly.graph_objects.Figure
        A plotly figure objects which can be displayed via the show() function,
        see https://plotly.com/python-api-reference/generated/plotly.graph_objects.Figure.html
    """

    logger = logging.getLogger("OpfLogger")

    n = graph_dict["N"]
    m = graph_dict["M"]

    adj = {}
    for i in range(m):
        adj[i] = [graph_dict[i][0] - 1, graph_dict[i][1] - 1]

    if m != alldata["numbranches"]:
        raise ValueError(
            "Number of branches in graph = %d, whereas number of branches in problem = %d."
            % (m, alldata["numbranches"])
        )

    if n != alldata["numbuses"]:
        raise ValueError(
            "Number of buses in graph = %d, whereas number of buses in problem = %d."
            % (n, alldata["numbuses"])
        )

    (
        vertexx,
        vertexy,
        scanned_list_consolidated,
        scanned_degrees_consolidated,
        scanned_unique_ordered_pairs,
        scanned_num_unique,
    ) = scangraph(alldata, graph_dict)

    # endbus is a list of end buses for network branches as ordered in graph_dict
    # revendbus is the reverse index
    # scanned_list_consolidated, degrees_consolidated are dictionaries whose keys are the vertex pairs (each pair in sorted order) corresponding to
    # lines graph file file.  list_consolidated is the array of lines whose ends are that vertex pair, listing the order in which
    # the lines are found in the graph file.  degrees_consolidated is the number of such lines (for each sorted vertex pair found in the file

    loud = False
    local_reordered_width = {}
    local_reordered_color = {}

    for j in range(scanned_num_unique):
        scannedordpair = scanned_unique_ordered_pairs[j]
        degscanned = scanned_degrees_consolidated[scannedordpair]
        if scannedordpair not in myedge_degrees_consolidated.keys():
            logger.info(
                "Error. Ordered pair %d -> (%d,%d) not in consolidated list\n"
                % (j, scannedordpair[0], scannedordpair[1])
            )
            raise ValueError(
                "Error. Ordered pair %d -> (%d,%d) not in consolidated list\n"
                % (j, scannedordpair[0], scannedordpair[1])
            )
        degmine = myedge_degrees_consolidated[scannedordpair]
        if degscanned != degmine:
            logger.info(
                "Error. Ordered pair %d -> (%d,%d) has different degrees %d, %d in scanned vs consolidated lists\n"
                % (j, scannedordpair[0], scannedordpair[1], degscanned, degmine)
            )
            raise ValueError(
                "Error. Ordered pair %d -> (%d,%d) has different degrees %d, %d in scanned vs consolidated lists\n"
                % (j, scannedordpair[0], scannedordpair[1], degscanned, degmine)
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
                logger.info(
                    " Width %d edge (%d,%d) orderings %d -> %d color %s.\n"
                    % (width, scannedordpair[0], scannedordpair[1], j, jprime, color)
                )

    if (
        True
    ):  # alldata['coordsfilename'] != None:  #should change this #TODO-Dan change to what?
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

    pos = {}
    for j in range(n):
        pos[j] = [vertexx[j], vertexy[j]]

    gG = Grbgraph()

    for j in range(n):
        gG.addvertex(j)

    newdeg = {}
    for j in range(scanned_num_unique):
        scannedordpair = scanned_unique_ordered_pairs[j]
        newdeg[scannedordpair[0], scannedordpair[1]] = 0

    reordered_width = {}
    reordered_color = {}
    reordered_position = {}

    loud = False
    for j in range(m):
        small = min(adj[j][0], adj[j][1])
        large = max(adj[j][0], adj[j][1])
        # G.add_edge(small, large)
        addcode = gG.addedge(small, large)
        if addcode:
            logger.info(
                "Could not add edge (%d, %d) to graph object.\n" % (small, large)
            )
            raise ValueError(
                "Could not add edge (%d, %d) to graph object.\n" % (small, large)
            )
        fbus = small + 1
        tbus = large + 1
        pair = (fbus, tbus)
        deg = newdeg[pair]
        if ((fbus, tbus)) in scanned_list_consolidated.keys():
            if loud and local_reordered_width[pair][deg] > 0:
                logger.info(
                    " pair ind0 %d (%d, %d) local color %s width %d\n"
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
    logger.info("Creating visualization object.\n")

    PH = Plotlyhandler(
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

    logger.info("Rendering figure.\n")

    xgap = np.max(vertexx) - np.min(vertexx)
    ygap = np.max(vertexy) - np.min(vertexy)

    myheight = 800
    mywidth = int(math.ceil(xgap / ygap * myheight))

    fig = PH.create_figure(height=myheight, width=mywidth, showlabel=False)
    logger.info("Created figure.\n")
    return fig


def scangraph(alldata, graph_dict):
    """
    Plots an OPF network out of user given solution/violation data which has been
    translated into necessary dictionaries in function grbgraphical.py:grbgraphical.

    Parameters
    ----------
    alldata : dictionary
        Main dictionary holding all necessary data
    graph_dict : dictionary
        Dictionary holding graph data

    Returns
    -------
    list
        List holding x coordinate values for each node
    list
        List holding y coordinate values for each node
    dictionary
        Dictionary holding possible multi-edges
    dictionary
        Dictionary holding the degree of each edge
    dictionary
        Dictionary holding each edge as a unique ordered pair
    int
        Number of unique edges

    """

    logger = logging.getLogger("OpfLogger")
    logger.info("Scanning graph.")
    logger.info("Using given lat, lon coordinates.\n")
    N = alldata["numbuses"]
    buses = alldata["buses"]
    nodex = np.zeros(N)
    nodey = np.zeros(N)
    for node in range(N):
        nodex[node] = buses[node + 1].lon
        nodey[node] = buses[node + 1].lat

    N = graph_dict["N"]
    M = graph_dict["M"]
    buses = alldata["buses"]

    trueM = 0
    endbus = {}
    revendbus = {}
    scanned_list_consolidated = {}
    scanned_degrees_consolidated = {}
    scanned_unique_ordered_pairs = {}
    scanned_num_unique = 0
    loud = False
    for i in range(M):
        edge = graph_dict[i]
        busfrom = edge[0]
        busto = edge[1]
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
        nodex,
        nodey,
        scanned_list_consolidated,
        scanned_degrees_consolidated,
        scanned_unique_ordered_pairs,
        scanned_num_unique,
    )


def grbgetgraphattr(alldata, value):
    # looks through threshold list to find first match

    color = "Black"
    size = 0
    caughtit = False
    numfeatures = 0  # alldata['graphical']['numfeatures']

    if numfeatures > 0:
        sizeval = alldata["graphical"]["sizeval"]
        colorstring = alldata["graphical"]["colorstring"]
        thresh = alldata["graphical"]["thresh"]

        size = 0  # sizeval[numfeatures-1]
        color = "Black"  # colorstring[numfeatures-1]

        for f in range(numfeatures):
            feature = numfeatures - f - 1
            # print(feature, 'value', value, thresh[feature], numfeatures)
            if value >= thresh[feature]:
                size = sizeval[feature]
                color = colorstring[feature]
                caughtit = True
                break

    else:  # if numfeatures == 0: #default case
        if value > 500:
            size = 20
            color = "orange"
            caughtit = True
        elif value > 150:
            size = 15
            color = "red"
            caughtit = True
        elif value > 100:
            size = 10
            color = "purple"
            caughtit = True

    return caughtit, size, color
