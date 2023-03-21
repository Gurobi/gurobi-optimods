import sys
import math
from log import danoLogger
from gurobipy import *
import numpy as np
from myutils import breakexit
import grbcasereader


def setconstants(log, all_data):
    log.joint("setting constants\n")
    epsilon_3 = 1e-3
    epsilon_2 = 1e-2

    all_data["epsilon_3"] = epsilon_3
    all_data["epsilon_2"] = epsilon_2

    log.joint("  epsilon_3: " + str(epsilon_3) + "\n")
    log.joint("  epsilon_2: " + str(epsilon_2) + "\n")
