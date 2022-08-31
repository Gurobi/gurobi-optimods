Maximum Weighted Matching
=========================

Something about specific graph algorithms vs mathprog formulations?

Problem Specification
---------------------

You are given a weighted graph G containing n vertices and m edges. Find the maximum weighted matching, i.e. select s subset of edges such that that no selected edge shares a vertex with any other selected edge, and the sum of edge weights for the selected set is maximised.

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
            \max \quad        & \sum_{(i, j) \in E} w_{ij} x_{ij} \\
            \mbox{s.t.} \quad & \sum_{j \in N(i)} x_{ij} \le 1 & \forall i \in V \\
                              & x_{ij} \in \lbrace 0, 1 \rbrace & \forall (i, j) \in E. \\
            \end{alignat}

        Note that for the general (non-bipartite) case, this model must be formulated and solved as a MIP, as there is no guarantee that simplex will return an integer solution for the relaxation.

|

Code
----

Tab between the code required to run the model from the store vs how to implement directly in gurobipy. If you use nupstup, all the gurobi internals are handled for you; users interact with the 'solver' by passing dataframes to a given spec and receiving a dataframe as output. If you instead peek under the hood and use gurobipy, you have more options to extend the model with additional constraints and data.

.. tabs::
    .. tab:: nupstup function

        .. literalinclude:: ../../../examples/weighted_matching/nupstup.py

    .. tab:: gurobipy model

        .. literalinclude:: ../../../examples/weighted_matching/gurobipy.py


Both codes construct the same model and give the same result. The model is solved as a LP/MIP/QP/etc by Gurobi.

.. collapse:: View Gurobi logs

    .. code-block:: text

        Gurobi Optimizer version 9.5.2 build v9.5.2rc0 (mac64[x86])
        Thread count: 4 physical cores, 8 logical processors, using up to 8 threads
        Optimize a model with 6 rows, 6 columns and 12 nonzeros
        Model fingerprint: 0xe6a49925
        Variable types: 0 continuous, 6 integer (6 binary)
        Coefficient statistics:
        Matrix range     [1e+00, 1e+00]
        Objective range  [1e+00, 1e+00]
        Bounds range     [1e+00, 1e+00]
        RHS range        [1e+00, 1e+00]
        Found heuristic solution: objective 2.2000000
        Presolve removed 6 rows and 6 columns
        Presolve time: 0.01s
        Presolve: All rows and columns removed

        Explored 0 nodes (0 simplex iterations) in 0.01 seconds (0.00 work units)
        Thread count was 1 (of 8 available processors)

        Solution count 2: 2.4 2.2 

        Optimal solution found (tolerance 1.00e-04)
        Best objective 2.400000000000e+00, best bound 2.400000000000e+00, gap 0.0000%

|

Solution
--------

Show the solution. Use doctests if possible (i.e. the solution must be stable enough). Otherwise, just display it somehow.

.. testcode:: nup
    :hide:

    from examples.weighted_matching.nupstup import matching

.. testoutput:: nup
    :hide:

    Gurobi Optimizer version ...
    ...

.. doctest:: nup
    :options: +NORMALIZE_WHITESPACE

    >>> matching
    <4x6 sparse matrix of type '<class 'numpy.float64'>'
        with 2 stored elements in COOrdinate format>
