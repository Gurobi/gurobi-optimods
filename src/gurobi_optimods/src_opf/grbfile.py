import sys
import os
import time
import math
import csv
import logging
import numpy as np


def initialize_data_dict(logfile):

    alldata = {}
    alldata["LP"] = {}
    alldata["MIP"] = {}
    alldata["logfile"] = logfile

    return alldata


def read_configfile(alldata, filename, casefile=""):
    """Function to read configurations for OPF solve from config file"""

    logging.info("Reading configuration file %s." % filename)
    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    # Currently Hard-coded
    alldata["usequadcostvar"] = False

    # TODO revisit default settings
    # Default values for options
    if casefile == "":
        casefilename = "NONE"
    else:
        casefilename = casefile
    dictionary_input = False
    voltsfilename = gvfilename = "NONE"
    lpfilename = "grbopf.lp"  # TODO default should be empty
    strictcheckvoltagesolution = (
        usevoltsolution
    ) = fixcs = skipjabr = cutplane = dodc = usemipstart = False
    useactivelossineqs = False
    useconvexformulation = False
    usemaxdispersion = False
    usemaxphasediff = False
    use_ef = False
    substitute_nonconv = False
    dographics = False
    graphattrsfilename = None
    coordsfilename = None
    dopolar = False
    doslp_polar = False
    doac = False
    doiv = False
    ivtype = None
    dodc = False
    branchswitching_mip = False
    branchswitching_comp = False
    maxdispersion_deg = 0
    maxphasediff_deg = 360
    fixtolerance = 0
    linenum = 0

    # Read lines of configuration file and save options
    while linenum < len(lines):
        thisline = lines[linenum].split()

        if len(thisline) <= 0:
            linenum += 1
            continue

        if thisline[0][0] == "#":
            linenum += 1
            continue

        if thisline[0] == "casefilename":
            casefilename = thisline[1]

        elif thisline[0] == "dictionary_input":
            dictionary_input = True

        elif thisline[0] == "lpfilename":
            lpfilename = thisline[1]

        elif thisline[0] == "use_ef":
            use_ef = True

        elif thisline[0] == "substitute_nonconv":
            substitute_nonconv = True

        elif thisline[0] == "dopolar":
            dopolar = True

        elif thisline[0] == "doslp_polar":
            doslp_polar = True

        elif thisline[0] == "branchswitching_mip":
            branchswitching_mip = True

        elif thisline[0] == "branchswitching_comp":
            branchswitching_comp = True

        elif thisline[0] == "usemaxdispersion":
            usemaxdispersion = True
            maxdispersion_deg = float(thisline[1])

        elif thisline[0] == "usemaxphasediff":
            usemaxphasediff = True
            maxphasediff_deg = float(thisline[1])

        elif thisline[0] == "strictcheckvoltagesolution":
            strictcheckvoltagesolution = True
            voltsfilename = thisline[1]

        elif thisline[0] == "voltsfilename":
            if len(thisline) < 3:
                raise ValueError(
                    "Voltsfilename option requires filename and tolerance."
                )

            voltsfilename = thisline[1]
            fixtolerance = float(thisline[2])
            if usevoltsolution:
                loging.error("Cannot use both voltsfilename and voltsolution\n")
                raise ValueError(
                    "Illegal option combination. Cannot use both voltsfilename and voltsolution."
                )

        elif thisline[0] == "usevoltsoution":
            if voltsfilename != "NONE":
                logging.error("Cannot use both voltsfilename and voltsolution\n")
                raise ValueError(
                    "Illegal option combination. Cannot use both voltsfilename and voltsolution."
                )

            usevoltsolution = True  # forcing solution in casefile

        elif thisline[0] == "FIXCS":
            fixcs = True

        elif thisline[0] == "useconvexformulation":
            useconvexformulation = True
            doac = True

        elif thisline[0] == "doiv":
            doiv = True
            use_ef = True  # needed
            doac = False
            dodc = False

            ivtype = thisline[1]

        elif thisline[0] == "gvfileoutput":
            gvfilename = thisline[1]

        elif thisline[0] == "skipjabr":
            skipjabr = True

        elif thisline[0] == "useactivelossineqs":
            useactivelossineqs = True

        elif thisline[0] == "usemipstart":
            usemipstart = True

        elif thisline[0] == "cutplane":
            cutplane = True

        elif thisline[0] == "dodc":
            dodc = 1
            doac = 0
            doiv = 0
            dodcbasic = 0

        elif thisline[0] == "dodcbasic":
            dodc = 0
            doac = 0
            dodcbasic = 1

        elif thisline[0] == "doac":
            dodc = 0
            doac = 1
            doiv = 0
            dodcbasic = 0

        elif thisline[0] == "dographics":
            dographics = True
            if len(thisline) > 1:
                graphattrsfilename = thisline[1]

        elif thisline[0] == "coordsfilename":
            coordsfilename = thisline[1]
        elif thisline[0] == "END":
            break

        else:
            raise ValueError("Illegal option input %s." % thisline[0])

        linenum += 1

    logging.info("Settings:")
    for x in [
        ("casefilename", casefilename),
        ("dictionary_input", dictionary_input),
        ("lpfilename", lpfilename),
        ("strictcheckvoltagesolution", strictcheckvoltagesolution),
        ("voltsfilename", voltsfilename),
        ("usevoltsolution", usevoltsolution),
        ("FIXCS", fixcs),
        ("useconvexformulation", useconvexformulation),
        ("skipjabr", skipjabr),
        ("useactivelossineqs", useactivelossineqs),
        ("usemipstart", usemipstart),
        ("cutplane", cutplane),
        ("dodc", dodc),
        ("doiv", doiv),
        ("ivtype", ivtype),
        ("doac", doac),
        ("fixtolerance", fixtolerance),
        ("use_ef", use_ef),
        ("substitute_nonconv", substitute_nonconv),
        ("dopolar", dopolar),
        ("doslp_polar", doslp_polar),
        ("usemaxdispersion", usemaxdispersion),
        ("maxdispersion_deg", maxdispersion_deg),
        ("usemaxphasediff", usemaxphasediff),
        ("maxphasediff_deg", maxphasediff_deg),
        ("dographics", dographics),
        ("graphattrsfilename", graphattrsfilename),
        ("coordsfilename", coordsfilename),
        ("branchswitching_mip", branchswitching_mip),
        ("branchswitching_comp", branchswitching_comp),
    ]:
        alldata[x[0]] = x[1]
        logging.info("  {} {}".format(x[0], x[1]))

    logging.info("")


def gbread_graphattrs(alldata, filename):
    """Function to read graphical attributes"""

    log = alldata["log"]

    log.joint("Reading graphical attributes file %s\n" % filename)
    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    numfeatures = 0

    thresh = {}
    colorstring = {}
    sizeval = {}

    linenum = 0
    # Read lines of configuration file and save options
    while linenum < len(lines):
        thisline = lines[linenum].split()

        if len(thisline) <= 0:
            linenum += 1
            continue

        if thisline[0] == "END":
            break

        thresh[numfeatures] = int(thisline[0])
        colorstring[numfeatures] = thisline[1]
        sizeval[numfeatures] = int(thisline[2])
        numfeatures += 1

        linenum += 1

    log.joint("Read %d graphical features\n" % numfeatures)
    alldata["graphical"]["numfeatures"] = numfeatures
    alldata["graphical"]["thresh"] = thresh
    alldata["graphical"]["colorstring"] = colorstring
    alldata["graphical"]["sizeval"] = sizeval


def grbread_coords(alldata):
    log = alldata["log"]

    # reads csv file with bus coordinates

    numbuses = alldata["numbuses"]
    buses = alldata["buses"]
    IDtoCountmap = alldata["IDtoCountmap"]

    filename = alldata["coordsfilename"]
    log.joint("reading csv coordinates file " + filename + "\n")

    f = open(filename, "r")
    csvf = csv.reader(f)
    thelist = list(csvf)
    f.close()

    for line in range(1, len(thelist)):
        # print(thelist[line][0], float(thelist[line][3]), float(thelist[line][4]))
        buses[line].lat = float(thelist[line][3])
        buses[line].lon = float(thelist[line][4])


def grbreadvoltsfile(alldata):
    log = alldata["log"]

    numbuses = alldata["numbuses"]
    buses = alldata["buses"]
    IDtoCountmap = alldata["IDtoCountmap"]

    filename = alldata["voltsfilename"]
    log.joint("reading volts file " + filename + "\n")

    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    linenum = 0

    while linenum < len(lines):
        thisline = lines[linenum].split()
        if len(thisline) > 0:

            if thisline[0] == "END":
                break
            elif thisline[0] == "bus":
                busid = int(thisline[1])
                buscount = IDtoCountmap[busid]

                vm = float(thisline[3])
                vadeg = float(thisline[5])
                varad = vadeg * math.pi / 180.0
                # print(busid, buscount, vm, vadeg, varad)

                buses[buscount].inputvoltage = 1
                buses[buscount].inputV = vm
                buses[buscount].inputA_rad = varad
            else:
                raise ValueError(
                    "Illegal line %d in voltsfile %s." % (linenum + 1, filename)
                )

        linenum += 1

    log.joint("Read volts file.\n")
