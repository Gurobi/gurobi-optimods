import logging
import math
import numpy as np

from gurobi_optimods.opf.grbgraph import Grbgraph
from gurobi_optimods.opf.plotlyhandler import Plotlyhandler

logger = logging.getLogger(__name__)


def generate_solution_figure(alldata, solution):
    """
    Generates a :class: `plotly.graph_objects.Figure` out of a given solution.
    Saves necessary data from solution to alldata dictionary
    and calls the main plotting function.

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param solution: Dictionary holding an OPF solution
    :type solution: dict

    :return: A plotly figure objects which can be displayed via the show() function,
             see https://plotly.com/python-api-reference/generated/plotly.graph_objects.Figure.html
    :rtype: :class: `plotly.graph_objects.Figure`
    """

    logger.info(
        f"Generating solution figure with objective value {solution['f']:.3f}. Coordinates given."
    )
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
    textlist.append([f"OBJ: {solution['f']:10.2f}", "black"])
    if numzeros > 0:
        textlist.append([f"Lines off: {numzeros}", "black"])
    else:
        textlist.append(["No lines turned off", "black"])
    textlist.append(["Black bus: generation <= 75 and load < 50", "black"])
    textlist.append(["Blue bus: generation <= 75 and load > 50", "blue"])
    textlist.append(["Purple bus: generation >  75", "purple"])
    textlist.append(["Orange bus: generation >  150", "orange"])
    textlist.append(["Red bus: generation >  500", "red"])

    return grbgraphical(alldata, "branchswitching", textlist)


def generate_violations_figure(alldata, violations):
    """
    Generates a :class: `plotly.graph_objects.Figure` out of given violations.
    Saves necessary data from solution to alldata dictionary
    and calls the main plotting function.

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param solution: Dictionary holding an OPF solution
    :type solution: dict

    :return: A plotly figure objects which can be displayed via the show() function,
             see https://plotly.com/python-api-reference/generated/plotly.graph_objects.Figure.html
    :rtype: :class: `plotly.graph_objects.Figure`
    """

    logger.info(f"Generating violations figure.")

    alldata["violation"] = {}
    alldata["violation"]["Vmagviol"] = {}  # Vm entry
    alldata["violation"]["IPviol"] = {}  # Pviol
    alldata["violation"]["IQviol"] = {}  # Qviol
    alldata["violation"]["branchlimit"] = {}  # limit violation
    numbuses = alldata["numbuses"]
    numbranches = alldata["numbranches"]

    maxvmviol = 0
    maxPviol = 0
    maxQviol = 0
    for i in range(1, numbuses + 1):
        busvmviol = violations["bus"][i]["Vmviol"]
        databus = alldata["buses"][i]
        alldata["violation"]["Vmagviol"][databus] = busvmviol
        maxvmviol = max(busvmviol, maxvmviol)

        busPviol = violations["bus"][i]["Pviol"]
        alldata["violation"]["IPviol"][databus] = busPviol
        maxPviol = max(busPviol, maxPviol)

        busQviol = violations["bus"][i]["Qviol"]
        alldata["violation"]["IQviol"][databus] = busQviol
        maxQviol = max(busPviol, maxQviol)

    maxlimiviol = 0
    for i in range(1, numbranches + 1):
        branchlimitviol = violations["branch"][i]["limitviol"]
        databranch = alldata["branches"][i]
        alldata["violation"]["branchlimit"][databranch] = branchlimitviol
        maxlimiviol = max(branchlimitviol, maxlimiviol)

    textlist = []
    textlist.append([f"max voltage magnitude violation: {maxvmviol:10.2f}", "black"])
    textlist.append([f"max real power injection violation: {maxPviol:10.2f}", "black"])
    textlist.append(
        [f"max reactive power injection violation: {maxQviol:10.2f}", "black"]
    )
    textlist.append([f"max branch limit violation: {maxlimiviol:10.2f}", "black"])
    textlist.append(
        ["Red bus: Vm violation > 1e-3 or (re)active power violation > 1e-2", "red"]
    )
    textlist.append(["Red branch: limit violation > 1e-3", "red"])
    return grbgraphical(alldata, "violation", textlist)


def grbgraphical(alldata, plottype, textlist):
    """
    Plots a given OPF solution. Depending on the plottype,
    the edge and node colors are different and different information
    is printed

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param plottype: Type of plot to show. Currently, `branchswitching` or `violation`
    :type plottype: str
    :param textlist: List of additional strings printed on the plot
    :type textlist: list

    :return: A plotly figure objects which can be displayed via the show() function,
             see https://plotly.com/python-api-reference/generated/plotly.graph_objects.Figure.html
    :rtype: :class: `plotly.graph_objects.Figure`
    """

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
    for branch in branches.values():
        graph_dict[counter] = (branch.count_f, branch.count_t)
        counter += 1

    node_text = {}
    mynode_size = {}
    mynode_color = {}
    mynode_border_width = {}
    edge_text = {}
    myedge_width = {}
    myedge_color = {}

    # Default actions
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
        Vmagviol = alldata["violation"]["Vmagviol"]  # Vm entry
        IPviol = alldata["violation"]["IPviol"]  # Pviol
        IQviol = alldata["violation"]["IQviol"]  # Qviol
        branchlimitviol = alldata["violation"]["branchlimit"]  # limit violation

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
            branch = branches[j]
            edge_text[j] = ""
            if abs(branchlimitviol[branch]) > 1e-3:
                myedge_width[j] = 5
                myedge_color[j] = "red"

    elif plottype == "branchswitching":

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

            largegen, mynode_size[j - 1], mynode_color[j - 1] = grbgetgraphattr(
                alldata, sumPgen
            )

            Pload = alldata["baseMVA"] * bus.Pd
            node_text[j - 1] = "Bus %d   Gen %7.2f  Load %7.2f" % (
                bus.nodeID,
                sumPgen,
                Pload,
            )
            # Default value for buses with load > 50 and that are not generating much
            if Pload > 50 and largegen == False:
                mynode_size[j - 1] = 7
                mynode_color[j - 1] = "blue"

        zholder = alldata["MIP"]["zholder"]
        for j in range(1, 1 + numbranches):
            branch = branches[j]
            f = branch.f
            t = branch.t
            count_of_f = IDtoCountmap[f]
            count_of_t = IDtoCountmap[t]

            edge_text[j] = "Branch {}   f {}  t {}".format(j, f, t)

            if zholder[j - 1] < 0.5:  # turned off
                myedge_width[j] = 5
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
        else:
            myedge_degrees_consolidated[(small, large)] += 1
            myedge_list_consolidated[(small, large)].append(j)

    return graphplot(
        alldata,
        graph_dict,
        node_text,
        mynode_size,
        mynode_color,
        mynode_border_width,
        edge_text,
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
    myedge_text,
    myedge_width,
    myedge_color,
    myedge_ends,
    myedge_list_consolidated,
    myedge_degrees_consolidated,
    textlist,
):
    """
    Generates a :class: `plotly.graph_objects.Figure` object from an OPF network out of user
    given solution/violation data which has been translated into necessary dictionaries
    in function ``grbgraphical``.
    Uses lat, lon coordinates as in coordinates dictionary,
    together with the plotlyhandler library to create a plotly figure.
    The figure is then rendered in a browser window.

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param graph_dict: Dictionary holding graph data
    :type graph_dict: dict
    :param myvertex_text: Dictionary holding the texts when hovering over a node
    :type myvertex_text: dict
    :param myvertex_size: Dictionary holding the size of each node. The size depends on generated power
    :type myvertex_size: dict
    :param myvertex_color: Dictionary holding the color of each node. The color depends on generated power
    :type myvertex_color: dict
    :param myvertex_border_width: Dictionary holding the thickness of the border of each node
    :type myvertex_border_width: dict
    :param myedge_width: Dictionary holding the thickness of each network edge
    :type myedge_width: dict
    :param myedge_color: Dictionary holding the color of each network edge. The color depends on violation and whether the edge is turned on/off
    :type myedge_color: dict
    :param myedge_ends: Dictionary holding each edge in both directions, e.g., (1,2) and (2,1). We need it, because, e.g., flows could be different
    :type myedge_ends: dict
    :param myedge_list_consolidated: Dictionary holding possible multi-edges. There could be parallel edges and we want to render all of them
    :type myedge_list_consolidated: dict
    :param myedge_degrees_consolidated: Dictionary holding the degree of each edge
    :type myedge_degrees_consolidated : dict
    :param textlist: List with additional strings printed on the plot
    :type textlist: list

    :raises ValueError: Inconsistent graph <-> data property

    :return: A plotly figure objects which can be displayed via the ``show()`` function,
             see https://plotly.com/python-api-reference/generated/plotly.graph_objects.Figure.html
    :rtype: :class: `plotly.graph_objects.Figure`
    """

    n = graph_dict["N"]
    m = graph_dict["M"]

    adj = {}
    for i in range(m):
        adj[i] = [graph_dict[i][0] - 1, graph_dict[i][1] - 1]

    if m != alldata["numbranches"]:
        raise ValueError(
            f"Number of branches in graph = {m}, whereas number of branches in problem = {alldata['numbranches']}."
        )

    if n != alldata["numbuses"]:
        raise ValueError(
            f"Number of buses in graph = {n}, whereas number of buses in problem = {alldata['numbuses']}."
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

    local_reordered_width = {}
    local_reordered_color = {}
    local_reordered_text = {}

    for j in range(scanned_num_unique):
        scannedordpair = scanned_unique_ordered_pairs[j]
        degscanned = scanned_degrees_consolidated[scannedordpair]
        if scannedordpair not in myedge_degrees_consolidated.keys():
            raise ValueError(
                f"Error. Ordered pair {j} -> ({scannedordpair[0]},{scannedordpair[1]}) not in consolidated list\n"
            )
        degmine = myedge_degrees_consolidated[scannedordpair]
        if degscanned != degmine:
            raise ValueError(
                f"Error. Ordered pair {j} -> ({scannedordpair[0]},{scannedordpair[1]}) has different degrees {degscanned}, {degmine} in scanned vs consolidated lists\n"
            )

        local_reordered_color[scannedordpair] = {}
        local_reordered_width[scannedordpair] = {}
        local_reordered_text[scannedordpair] = {}

        for k in range(degmine):
            j = myedge_list_consolidated[scannedordpair][k]
            jprime = scanned_list_consolidated[scannedordpair][k]
            color = myedge_color[j]
            width = myedge_width[j]
            text = myedge_text[j]

            local_reordered_color[scannedordpair][k] = color
            local_reordered_width[scannedordpair][k] = width
            local_reordered_text[scannedordpair][k] = text

    vertexx -= np.min(vertexx)
    vertexy -= np.min(vertexy)

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

    # Construct network graph
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
    reordered_text = {}

    for j in range(m):
        small = min(adj[j][0], adj[j][1])
        large = max(adj[j][0], adj[j][1])
        # G.add_edge(small, large)
        error = gG.addedge(small, large)
        if error:
            raise ValueError(
                f"Could not add edge ({small}, {large}) to graph object.\n"
            )
        fbus = small + 1
        tbus = large + 1
        pair = (fbus, tbus)
        deg = newdeg[pair]
        if ((fbus, tbus)) in scanned_list_consolidated.keys():
            newdeg[fbus, tbus] += 1

        reordered_width[j] = local_reordered_width[pair][deg]
        reordered_color[j] = local_reordered_color[pair][deg]
        reordered_text[j] = local_reordered_text[pair][deg]
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
        edge_position=reordered_position,
        edge_text=myedge_text,
        vertex_text=myvertex_text,  # should be reordered, for both?
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
    Scans given graph data and translates it into internal lists and dictionaries

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param graph_dict: Dictionary holding graph data
    :type graph_dict: dict

    :return: List holding x coordinate values for each node
             List holding y coordinate values for each node
             Dictionary holding possible multi-edges
             Dictionary holding the degree of each edge
             Dictionary holding each edge as a unique ordered pair
             Number of unique edges
    :rtype: list, list, dict, dict, dict, int
    """

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
        else:
            scanned_degrees_consolidated[(small, large)] += 1
            scanned_list_consolidated[(small, large)].append(trueM)

        endbus[trueM] = (busfrom, busto)
        revendbus[(busfrom, busto)] = trueM + 1
        revendbus[(busto, busfrom)] = trueM + 1
        trueM += 1
    logger.info(f"After scanning, number of edges is {trueM}.\n")
    return (
        nodex,
        nodey,
        scanned_list_consolidated,
        scanned_degrees_consolidated,
        scanned_unique_ordered_pairs,
        scanned_num_unique,
    )


def grbgetgraphattr(alldata, value):
    """
    Looks through a possibly pre-set graphical attributes list
    and sets the corresponding values in ``alldata`` dictionary

    :param alldata: Main dictionary holding all necessary data
    :type alldata: dict
    :param value: Threshold value for graph attributes
    :type value: float

    :return: States whether the graph attributes have been set, i.e.,
             if the value is > 100
             Size of node
             Color of node
    :rtype: bool, float, float
    """

    color = "Black"
    size = 2
    valuesset = False
    numfeatures = 0  # alldata['graphical']['numfeatures']

    # Currently unused, may be added if we decide to allow users to provide
    # their own graph attributes
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
                valuesset = True
                break

    else:  # if numfeatures == 0: #default case
        if value > 500:
            size = 16
            color = "red"
            valuesset = True
        elif value > 150:
            size = 13
            color = "orange"
            valuesset = True
        elif value > 75:
            size = 8
            color = "purple"
            valuesset = True

    return valuesset, size, color
