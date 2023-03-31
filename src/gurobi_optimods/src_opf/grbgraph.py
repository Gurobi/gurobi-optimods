import logging


class grbGraph:
    def __init__(self):
        self.vertices = {}
        self.edges = {}
        self.n = 0  # number of vertices
        self.m = 0  # number of edges
        logging.info("Created grbGraph object.")

    def addvertex(self, i):
        self.vertices[self.n] = i
        self.n += 1

    def addedge(self, i, j):
        if i in self.vertices.values() and j in self.vertices.values():
            self.edges[self.m] = (i, j)
            self.m += 1
            return 0
        else:
            return 1

    def getmetrics(self):
        logging.info("Graph object has %d vertices %d edges." % (self.n, self.m))
