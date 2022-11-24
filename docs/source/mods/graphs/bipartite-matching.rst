Maximum Cardinality Bipartite Matching
======================================

Something about specific graph algorithms vs mathprog formulations?

Problem Specification
---------------------

You are given a bipartite graph G containing n vertices and m edges. Find the maximum matching, i.e. select as many edges as possible so that no selected edge shares a vertex with any other selected edge.

.. tabs::

    .. tab:: Domain-Specific Description

        A matching is a set of pairwise non-adjacent edges ...

    .. tab:: Mathematical Model

        For each edge :math:`(i,j)` in graph :math:`G = (V, E)`, we should select a subset. So define variables as follows

        .. math::

            x_{ij} = \begin{cases}
                1 & \text{if edge}\,(i,j)\,\text{is selected in the matching} \\
                0 & \text{otherwise.} \\
            \end{cases}

        and the mathematical model as (check undirected terminology/notation here)

        .. math::

            \begin{alignat}{2}
            \max \quad        & \sum_{(i, j) \in E} x_{ij} \\
            \mbox{s.t.} \quad & \sum_{j \in N(i)} x_{ij} \le 1 & \forall i \in V \\
                              & 0 \le x_{ij} \le 1 & \forall (i, j) \in E \\
            \end{alignat}

        Note that for the bipartite case, simplex is sufficient, we do not need to use binary variables, just bounded ones.

|

Code
----

Show the code required to run the model from the store vs how to implement directly in gurobipy. All the gurobi internals are handled for you; users interact with the 'solver' by passing dataframes to a given spec and receiving a dataframe as output.

.. literalinclude:: ../../../examples/bipartite_matching.py

Both codes construct the same model and give the same result. The model is solved as a LP/MIP/QP/etc by Gurobi.

.. collapse:: View Gurobi logs

    .. code-block:: text

        Gurobi Optimizer version 9.5.2 build v9.5.2rc0 (mac64[x86])
        Thread count: 4 physical cores, 8 logical processors, using up to 8 threads
        Optimize a model with 7 rows, 6 columns and 12 nonzeros
        Model fingerprint: 0x66da6fea
        Coefficient statistics:
        Matrix range     [1e+00, 1e+00]
        Objective range  [1e+00, 1e+00]
        Bounds range     [0e+00, 0e+00]
        RHS range        [1e+00, 1e+00]
        Presolve removed 7 rows and 6 columns
        Presolve time: 0.01s
        Presolve: All rows and columns removed
        Iteration    Objective       Primal Inf.    Dual Inf.      Time
            0    3.0000000e+00   0.000000e+00   0.000000e+00      0s

        Solved in 0 iterations and 0.01 seconds (0.00 work units)
        Optimal objective  3.000000000e+00

|

Solution
--------

Show the solution. Use doctests if possible (i.e. the solution must be stable enough). Otherwise, just display it somehow.

.. testcode:: bipartite_matching
    :hide:

    from examples.bipartite_matching import matching, G

.. testoutput:: bipartite_matching
    :hide:

    ...
    Optimal objective  3.000000000e+00

.. doctest:: bipartite_matching
    :options: +NORMALIZE_WHITESPACE

    >>> matching
    <8x8 sparse matrix of type '<class 'numpy.float64'>'
        with 3 stored elements in COOrdinate format>

.. doctest:: bipartite_matching
    :options: +NORMALIZE_WHITESPACE

    >>> import networkx as nx
    >>> import matplotlib.pyplot as plt
    >>> g = nx.from_scipy_sparse_array(G)
    >>> layout = nx.bipartite_layout(g, [0, 1, 2, 3, 4])
    >>> fig, (ax1, ax2) = plt.subplots(1, 2)
    >>> nx.draw(g, layout, ax=ax1)
    >>> g = nx.from_scipy_sparse_array(matching)
    >>> nx.draw(g, layout, ax=ax2)

.. image:: bipartite-result.png
  :width: 600
  :alt: Bipartite matching result
