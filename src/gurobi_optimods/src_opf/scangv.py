import sys
import os
import numpy as np
import re

if __name__ == "__main__":
    filename = sys.argv[1]
    print(filename)

    try:
        f = open(filename, "r")
        lines = f.readlines()
        f.close()
    except:
        sys.exit("failure")

    N = len(lines)
    # thestrings = np.array(lines[linenum].split() for linenum in range(N))

    criterion = np.array([len(lines[linenum].split()) for linenum in range(N)])
    selection = np.where(criterion == 2)
    N = len(selection[0])
    node = np.full(N, -1)
    nodex = np.zeros(N)
    nodey = np.zeros(N)
    trueN = 0

    for k in range(N):
        i = selection[0][k]
        foo = lines[i].split()
        # print(k, 'i',i, 'split', foo[0], foo[1])
        if len(foo[1]) > 3 and foo[1][:4] == "[pos":
            # print(k, 'i',i, 'split', foo[0], foo[1])
            comma = foo[1].rfind(",")
            rightbr = foo[1].rfind("]")
            # print(comma, rightbr)
            node = int(foo[0]) - 1  # base 0
            nodex[node] = float(foo[1][6:comma])
            nodey[node] = float(foo[1][comma + 1 : rightbr - 1])
            # print(node, float(foo[1][6:comma]), float(foo[1][comma+1:rightbr-1]))
            trueN += 1
    for i in range(N):
        print(i, nodex[i], nodey[i])
