Network Flow
============

Network flow problems appear are a generalised problem type defined on a graph.
A common example of these would be the: min-cost flow, the max-flow/min-cut, the
shortest path problem.

Problem Specification
---------------------

Give a brief overview of the problem being solved.

.. tabs::

    .. tab:: Graph Theory

        Give a definition of the problem in the language of the domain expert.
        For a given graph :math:`G` with set of vertices :math:`V` and edges :math:`E`. Depending on the problem, each edge may have differnt parameters associated to it for instance, for a given :math:`(i,j)\in E`, we may have costs :math:`c_{ij}`, capacities :math:`B_{ij}` and each node :math:`i` may have a demand :math:`d_i`.

        - Min-cost flow: minimize the cost subject to edge capacities and fulfil
          node demand.
        - Max-flow/min-cut: maximize flow through the graph subject to edge
          capacities.
        - Shortest path: find a path from source to sink that has the minimal
          cost.

    .. tab:: Optimization Model

        For a given graph :math:`G` with set of vertices :math:`V` and edges :math:`E`. Depending on the problem, each edge may have differnt parameters associated to it for instance, for a given :math:`(i,j)\in E`, we may have costs :math:`c_{ij}`, capacities :math:`B_{ij}` and each node :math:`i` may have a demand :math:`d_i`.

        We need a set of continuous variables :math:`x_{ij}` to denote the
        amount of flow going through an edge.

        .. math::

            \begin{alignat}{2}
              \min \quad        & \sum_{(i, j) \in E} c_{ij} x_{ij} \\
              \mbox{s.t.} \quad & \sum_{j \in \delta^+(i)} x_{ij} - \sum_{j \in \delta^-(i)} x_{ji} = d_i & \forall i \in V' \\
                                & 0 \leq x_{ij} \le B_{ij} & \forall (i, j) \in E \\
            \end{alignat}

        Where :math:`\delta^+(\cdot)` (:math:`\delta^-(\cdot)`) denotes the outgoing (incoming) neighours.



.. tabs::

    .. tab:: ``load_network_flow_example_data``

        Let us load a sample `networkx` graph.

        .. doctest:: mod
            :options: +NORMALIZE_WHITESPACE

            >>> from gurobi_optimods import datasets
            >>> G, s, t = datasets.load_network_flow_example_data()
            >>> for n in G.nodes(data=True):
            ...     print(n)
            ...
            (0, {'demand': 20})
            (1, {'demand': 0})
            (2, {'demand': 0})
            (3, {'demand': -5})
            (4, {'demand': -15})
            >>> for e in G.edges(data=True):
            ...     print(e)
            ...
            (0, 1, {'capacity': 15, 'cost': 4})
            (0, 2, {'capacity': 8, 'cost': 4})
            (1, 3, {'capacity': 4, 'cost': 2})
            (1, 2, {'capacity': 20, 'cost': 2})
            (1, 4, {'capacity': 10, 'cost': 6})
            (2, 3, {'capacity': 15, 'cost': 1})
            (2, 4, {'capacity': 5, 'cost': 3})
            (3, 4, {'capacity': 20, 'cost': 2})
            (4, 2, {'capacity': 4, 'cost': 3})
            >>> s, t
            (0, 4)

    .. tab:: Graph

      TODO: add plot

|

Code
----

Self contained code example to run the mod from an example dataset. Example
datasets should bd included in the ``gurobi_optimods.datasets`` module for
easy access by users.

.. testcode:: mod

    from gurobi_optimods.datasets import load_network_flow_example_data
    from gurobi_optimods.network_flow import solve_min_cost


    G, source, sink = load_network_flow_example_data()
    sol, cost = shortest_path(G, source, sink)

..  A snippet of the Gurobi log output here won't show in the rendered page,
    but serves as a doctest to make sure the code example runs. The ... lines
    are meaningful here, they will match anything in the output test.

.. testoutput:: mod
    :hide:
    ...
    Solved in 2 iterations and 0.00 seconds (0.00 work units)
    Optimal objective  7.000000000e+00
    ...

The model is solved as an LP by Gurobi.

.. collapse:: View Gurobi Logs

    .. code-block:: text

      Gurobi Optimizer version 10.0.1 build v10.0.1rc0 (mac64[arm])

      CPU model: Apple M1
      Thread count: 8 physical cores, 8 logical processors, using up to 8 threads

      WLS license - registered to david.torres-sanchez@gurobi.com
      Optimize a model with 5 rows, 9 columns and 18 nonzeros
      Model fingerprint: 0xce040fba
      Coefficient statistics:
        Matrix range     [1e+00, 1e+00]
        Objective range  [1e+00, 6e+00]
        Bounds range     [1e+00, 1e+00]
        RHS range        [1e+00, 1e+00]
      Presolve removed 1 rows and 1 columns
      Presolve time: 0.00s
      Presolved: 4 rows, 8 columns, 16 nonzeros

      Iteration    Objective       Primal Inf.    Dual Inf.      Time
             0    4.0000000e+00   2.000000e+00   0.000000e+00      0s
             2    7.0000000e+00   0.000000e+00   0.000000e+00      0s

      Solved in 2 iterations and 0.00 seconds (0.00 work units)
      Optimal objective  7.000000000e+00

|

Solution
--------

Show the solution. One way is to use doctests to display simple shell outputs
(see the workforce example). This can be done simply by pasting outputs
directly from a python shell. Another option is to include and display figures
(see the graph matching examples).

.. doctest:: mod
    :options: +NORMALIZE_WHITESPACE

    >>>
