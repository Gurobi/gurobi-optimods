import sys
import os
import numpy as np
import re

def scangv(filename):
    print('Scanning gv file',filename)

    try:
        f     = open(filename, "r")
        lines = f.readlines()
        f.close()
    except:
        return 0, None, None

    N = len(lines)
    #thestrings = np.array(lines[linenum].split() for linenum in range(N))

    criterion = np.array([ len(lines[linenum].split()) for linenum in range(N) ])
    selection = np.where(criterion == 2)
    N         = len(selection[0])
    node      = np.full(N,-1)
    nodex     = np.zeros(N)
    nodey     = np.zeros(N)
    trueN     = 0

    for k in range(N):
        i   = selection[0][k]
        foo = lines[i].split()
        #print(k, 'i',i, 'split', foo[0], foo[1])
        if len(foo[1]) > 3 and foo[1][:4]== '[pos':
            #print(k, 'i',i, 'split', foo[0], foo[1])
            #print(foo)
            comma   = foo[1].rfind(',')
            rightbr = foo[1].rfind(']')
            #print(comma, rightbr)
            node        = int(foo[0])-1 #base 0
            nodex[node] = float(foo[1][6:comma])
            nodey[node] = float(foo[1][comma+1:rightbr-1])
            #print(node, float(foo[1][6:comma]), float(foo[1][comma+1:rightbr-1]))
            trueN += 1

    thisM = 0
    selection = np.where(criterion > 3)
    N = len(selection[0])
    endbus = {}
    for k in range(N):
        i   = selection[0][k]
        foo = lines[i].split()
        if foo[1][0] == '-':
            #print(len(foo),foo)
            busfrom = int(foo[0])
            busto = int(foo[2])
            #print(busfrom, busto)
            endbus[thisM] = (busfrom, busto)
            thisM += 1
    

    return trueN, N, nodex, nodey, thisM, endbus
