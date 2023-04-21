Maximum Bipartite Matching
==========================

Something about specific graph algorithms vs mathprog formulations?

Problem Specification
---------------------

Consider a bipartite graph :math:`G(U, V, E)`, where :math:`U` and :math:`V`
are disjoint vertex sets, and the edge set :math:`E \subseteq U \times V`
joins only between, not within, the sets. A matching on this graph is any
subset of edges such that no vertex is incident to more than one edge. A
maximum matching is the largest possible matching on :math:`G`.

The bipartite matching problem can be reduced to a maximum flow problem

You are given a bipartite graph :math:`G` containing :math:`n` vertices and
:math:`m` edges. Find the maximum matching, i.e. select as many edges as
possible so that no selected edge shares a vertex with any other selected edge.

.. tabs::

    .. tab:: Domain-Specific Description

        A matching is a set of pairwise non-adjacent edges ...

    .. tab:: Optimization Model

        Use a network max-flow model ...

        Note that for the bipartite case, simplex is sufficient, we do not
        need to use binary variables, just :math:`[0, 1]` bounds. Gurobi uses
        a network simplex algorithm to solve such models.

|

Code
----

.. testcode:: bipartite_matching

    import numpy as np
    import scipy.sparse as sp

    from gurobi_optimods.matching import maximum_bipartite_matching

    # Create a simple bipartite graph as a sparse matrix
    nodes1 = np.array([0, 1, 2, 3, 4])
    nodes2 = np.array([5, 6, 7])
    row = [0, 3, 4, 0, 1, 3]
    col = [7, 5, 5, 6, 6, 7]
    data = [1, 1, 1, 1, 1, 1]
    adjacency = sp.coo_array((data, (row, col)), shape=(8, 8))

    # Compute the maximum matching
    matching = maximum_bipartite_matching(adjacency, nodes1, nodes2)

.. testoutput:: bipartite_matching
    :hide:

    ...
    Optimal objective -3.000000000e+00

The `maximum_bipartite_matching` function formulates a linear program for the
the network flow model corresponding to the given bipartite graph. Gurobi
solves this model using a network primal simplex algorithm.

.. collapse:: View Gurobi logs

    .. code-block:: text

        Solving maximum matching n1=5 n2=3 |E|=6
        Maximum matching formulated as min-cost flow with 10 nodes and 15 arcs
        Restricted license - for non-production use only - expires 2024-10-28
        Gurobi Optimizer version 10.0.1 build v10.0.1rc0 (mac64[x86])

        CPU model: Intel(R) Core(TM) i5-1038NG7 CPU @ 2.00GHz
        Thread count: 4 physical cores, 8 logical processors, using up to 8 threads

        Optimize a model with 10 rows, 15 columns and 30 nonzeros
        Model fingerprint: 0xb08809c2
        Coefficient statistics:
          Matrix range     [1e+00, 1e+00]
          Objective range  [1e+00, 1e+00]
          Bounds range     [1e+00, 1e+00]
          RHS range        [0e+00, 0e+00]
        Presolve removed 4 rows and 4 columns
        Presolve time: 0.00s
        Presolved: 6 rows, 11 columns, 22 nonzeros

        Iteration    Objective       Primal Inf.    Dual Inf.      Time
               0   -3.0000000e+00   1.000000e+00   0.000000e+00      0s
               1   -3.0000000e+00   0.000000e+00   0.000000e+00      0s

        Solved in 1 iterations and 0.00 seconds (0.00 work units)
        Optimal objective -3.000000000e+00
        Done: max bipartite matching has 3 edges

|

Solution
--------

The maximum matching is returned as a subgraph of the original bipartite
graph, as a ``scipy.sparse`` array. Inspecting the result, it is clear that
this is a maximum matching, since no two edges share a node in common, and
all nodes in the second set are incident to an edge in the matching.

.. doctest:: bipartite_matching
    :options: +NORMALIZE_WHITESPACE

    >>> print(sp.triu(matching))
      (0, 7)        1.0
      (1, 6)        1.0
      (3, 5)        1.0

We can also inspect the result by plotting the graph and the edges selected
in the matching using networkx.

.. doctest:: bipartite_matching
    :options: +NORMALIZE_WHITESPACE

    >>> import networkx as nx
    >>> import matplotlib.pyplot as plt
    >>>
    >>> fig, (ax1, ax2) = plt.subplots(1, 2)
    >>> g = nx.from_scipy_sparse_array(adjacency)
    >>> layout = nx.bipartite_layout(g, nodes1)
    >>> nx.draw(g, layout, ax=ax1)
    >>> g = nx.from_scipy_sparse_array(matching)
    >>> nx.draw(g, layout, ax=ax2)

.. image:: figures/bipartite-result.png
  :width: 600
  :alt: Bipartite matching result
