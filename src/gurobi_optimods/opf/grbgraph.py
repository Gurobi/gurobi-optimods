import logging


class Grbgraph:
    """
    Simple network graph class to hold necessary edge and vertex information.
    """

    def __init__(self):
        self.vertices = {}
        self.edges = {}
        self.n = 0  # number of vertices
        self.m = 0  # number of edges
        logger = logging.getLogger("OpfLogger")
        logger.info("Created Grbgraph object.")

    def addvertex(self, i):
        """
        Adds a vertex to the Graph

        :param i: vertex index
        :type i: int
        """
        self.vertices[self.n] = i
        self.n += 1

    def addedge(self, i, j):
        """
        Adds an edge to the Graph

        :param i: vertex index of first node
        :type i: int
        :param j: vertex index of second node
        :type j: int

        :return: 0 if edge could be added, 1 otherwise
        :rtype: int
        """
        if i in self.vertices.values() and j in self.vertices.values():
            self.edges[self.m] = (i, j)
            self.m += 1
            return 0
        else:
            return 1

    def getmetrics(self):
        """
        Prints graph statistics
        """
        logger = logging.getLogger("OpfLogger")
        logger.info(f"Graph object has {self.n} vertices {self.m} edges.")
