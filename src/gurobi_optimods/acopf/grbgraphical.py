import sys
import math
#from log import Logger
from gurobipy import *
import numpy as np
from myutils import breakexit
import time
from graph4 import *

epsilon4 = 1e-4

def grbgraphical(alldata):
  log = alldata['log']
  log.joint("graphical layout, 2\n")

  numbuses = alldata['numbuses']
  numbranches = alldata['numbranches']

  gvfilename = 'grbgraphical.gv'
  txtfilename = 'newgraph.txt'

  try:
    f = open(gvfilename, "w")
    log.joint("writing to gv file " + gvfilename + "\n")
    g = open(txtfilename, "w")
  except:
    log.stateandquit("cannot open file " + gvfilename)
  f.write("graph {\n")
  f.write('node [color=black, height=0, label=\"\\N\", shape=point, width=0];\n')
  g.write('N '+str(numbuses) + ' M ' + str(numbranches) + '\n')
  for bus in alldata['buses'].values():
    #f.write("     " + str(bus.nodeID)+";\n")
    f.write("     " + str(bus.count)+";\n")
  for branch in alldata['branches'].values():
    f.write("     " + str(branch.id_f)+" -- " + str(branch.id_t)+";\n")
    g.write(' ' + str(branch.id_f)+ ' ' + str(branch.id_t)+ '\n')

  f.write("}\n")
  f.close()
  g.write('END\n')
  g.close()

  scale = 10
  firstgvfile = 'first.gv'
  sfdpcommand = 'sfdp -Goverlap_scaling='+str(scale)+' -o '+firstgvfile + ' ' + gvfilename
  log.joint('sfdp command: ' + sfdpcommand + '\n')
  system(sfdpcommand)

  '''
  jpgfile = 'second.jpg'
  neatocommand = 'neato -Tjpeg -n -o '+ jpgfile + ' ' + firstgvfile
  log.joint('neato command: ' + neatocommand + '\n')
  system(neatocommand)
  '''
  
  breakexit('graph,1')

  graphplot(txtfilename, firstgvfile)

  breakexit('graph,2')
