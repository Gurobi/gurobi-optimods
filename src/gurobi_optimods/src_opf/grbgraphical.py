import sys
import math
import time
import logging
from .graph4 import *
from gurobipy import *
from .utils import break_exit


def plot_solution(alldata, solution, objval):
    """Plot feasible solution"""

    logger = logging.getLogger("OpfLogger")
    logger.info("Plotting solution with value %.3e. Coordinates given." % objval)
    x = list(solution.values())
    numbuses = alldata["numbuses"]
    buses = alldata["buses"]
    numbranches = alldata["numbranches"]
    branches = alldata["branches"]
    gens = alldata["gens"]
    IDtoCountmap = alldata["IDtoCountmap"]

    # thetavar = alldata["LP"]["thetavar"]
    # Pvar_f = alldata["LP"]["Pvar_f"]
    # twinPvar_f = alldata["LP"]["twinPvar_f"]

    zholder = np.zeros(numbranches)
    alldata["MIP"]["zholder"] = zholder
    gholder = np.zeros(alldata["numgens"])
    alldata["MIP"]["gholder"] = gholder
    zholder = alldata["MIP"]["zholder"]
    gholder = alldata["MIP"]["gholder"]

    numzeros = 0
    for j in range(1, 1 + numbranches):
        branch = branches[j]
        zholder[j - 1] = x[branch.switchvarind]
        if x[branch.switchvarind] < 0.5:
            numzeros += 1

    for j1 in range(1, 1 + alldata["numgens"]):
        gen = alldata["gens"][j1]
        gholder[j1 - 1] = (
            alldata["baseMVA"] * x[gen.Pvarind]
        )  # print('gen',gen.count,x[gen.Pvarind])
        # break_exit('printed')  #functionality to be added # TODO-Dan What functionality is that?
    textlist = []
    textlist.append("OBJ: %10.2f" % (objval))
    textlist.append("Lines off: %d" % (numzeros))
    grbgraphical(alldata, "branchswitching", textlist)


def grbgraphical(alldata, plottype, textlist):
    logger = logging.getLogger("OpfLogger")
    buses = alldata["buses"]
    numbuses = alldata["numbuses"]
    branches = alldata["branches"]
    numbranches = alldata["numbranches"]
    gens = alldata["gens"]
    IDtoCountmap = alldata["IDtoCountmap"]
    txtfilename = "newgraph.txt"
    logger.info("Graphical layout, coordinates given.\n")

    g = open(txtfilename, "w")
    logger.info("Writing to txt file %s\n" % txtfilename)

    g.write("N " + str(numbuses) + " M " + str(numbranches) + "\n")
    for branch in alldata["branches"].values():
        g.write(" " + str(branch.count_f) + " " + str(branch.count_t) + "\n")
    g.write("END\n")
    g.close()

    # break_exit('graph,1')
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
        branch = alldata["branches"][j]
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
            node_text[j - 1] = "Bus %d   Gen %7.2f  Load %7.2f" % (j, sumPgen, Pload)

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

    graphplot(
        alldata,
        txtfilename,
        node_text,
        mynode_size,
        mynode_color,
        mynode_border_width,
        myedge_width,
        myedge_color,
        myedge_ends,
        myedge_list_consolidated,
        myedge_degrees_consolidated,
        numbranches,
        textlist,
    )
    # break_exit('graph,2')


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

    if False:  # caughtit:
        print(value, caughtit, size, color)
        break_exit("caughtit")

    return caughtit, size, color
