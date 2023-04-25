Minimum Cost Flow
=================

Minimum cost flow problems are defined on a graph or network where the goal is
to send an amount of flow. Other graph problems can be modelled using this
framework (with suitable alterations).

Problem Specification
---------------------


.. tabs::

    .. tab:: Graph Theory

        For a given graph :math:`G` with set of vertices :math:`V` and edges
        :math:`E`. For a given :math:`(i,j)\in E`, we have:

        - cost: :math:`c_{ij}\in \mathbb{R}`;
        - and capacity: :math:`B_{ij}\in\mathbb{R}`.

        Also, each vertex :math:`i\in V` has a demand :math:`d_i\in\mathbb{R}`.

        The problem can be stated as finding a capacity feasible flow such that
        the total cost is minimised and the demand at each vertex is met.

    .. tab:: Optimization Model

        For a given graph :math:`G` with set of vertices :math:`V` and edges
        :math:`E`. For a given :math:`(i,j)\in E`, we have:

        - cost: :math:`c_{ij}\in \mathbb{R}`;
        - and capacity: :math:`B_{ij}\in\mathbb{R}`.

        Also, each vertex :math:`i\in V` has a demand :math:`d_i\in\mathbb{R}`.

        Let us define a set of continuous variables :math:`x_{ij}` to represent
        the amount of flow going through an edge.


        The mathematical formulation can be stated as follows:

        .. math::

            \begin{alignat}{2}
              \min \quad        & \sum_{(i, j) \in E} c_{ij} x_{ij} \\
              \mbox{s.t.} \quad & \sum_{j \in \delta^+(i)} x_{ij} - \sum_{j \in \delta^-(i)} x_{ji} = d_i & \forall i \in V \\
                                & 0 \leq x_{ij} \le B_{ij} & \forall (i, j) \in E \\
            \end{alignat}

        Where :math:`\delta^+(\cdot)` (:math:`\delta^-(\cdot)`) denotes the outgoing (incoming) neighours.

The input data for this mod includes interfaces for the following:

* Pandas DataFrame with the format described below.
* A networkx graph with the edge/vertex attributes previously mentioned.
* Scipy sparse matrices in CSR format containing the graph structure, and the capacity and cost matrix. Additionally a numpy array representing the demands of the vertices.

|

Code
----

.. testcode:: network_flow

    from gurobi_optimods.datasets import load_min_cost_flow
    from gurobi_optimods.min_cost_flow import min_cost_flow


    edge_data, node_data = load_min_cost_flow()
    cost, sol = min_cost_flow(edge_data, node_data)

.. testoutput:: network_flow
    :hide:

    ...
    Optimal objective  1.500000000e+02

The model is solved as an LP by Gurobi.

.. collapse:: View Gurobi Logs

    .. code-block:: text

      Solving min-cost flow with 5 nodes and 9 edges
      Gurobi Optimizer version 10.0.1 build v10.0.1rc0 (mac64[arm])

      CPU model: Apple M1
      Thread count: 8 physical cores, 8 logical processors, using up to 8 threads

      WLS license - registered to david.torres-sanchez@gurobi.com
      Optimize a model with 5 rows, 9 columns and 18 nonzeros
      Model fingerprint: 0x516e531a
      Coefficient statistics:
        Matrix range     [1e+00, 1e+00]
        Objective range  [1e+00, 6e+00]
        Bounds range     [4e+00, 2e+01]
        RHS range        [5e+00, 2e+01]
      Presolve removed 2 rows and 2 columns
      Presolve time: 0.00s
      Presolved: 3 rows, 7 columns, 14 nonzeros

      Iteration    Objective       Primal Inf.    Dual Inf.      Time
             0    7.2994000e+01   5.000400e+01   0.000000e+00      0s
             2    1.5000000e+02   0.000000e+00   0.000000e+00      0s

      Solved in 2 iterations and 0.00 seconds (0.00 work units)
      Optimal objective  1.500000000e+02

|

Solution
--------

.. tabs::

    .. tab:: ``load_min_cost_flow``

        Let us load a sample graph in a pandas DataFrame format.

        .. doctest:: mod
            :options: +NORMALIZE_WHITESPACE

            >>> from gurobi_optimods import datasets
            >>> edge_data, node_data = datasets.load_min_cost_flow2()
            >>> edge_data
                           capacity  cost
            source target
            0      1              2     9
                   2              2     7
            1      3              1     1
            2      3              1    10
                   4              2     6
            3      5              2     1
            4      5              2     1
            >>> node_data
               posx  posy  demand original_label
            0     0   0.0      -2              s
            1     1   0.5       0              1
            2     1  -0.5      -1              2
            3     2   0.5       1              3
            4     2  -0.5       0              4
            5     3   0.0       2              t

    .. tab:: ``networkx`` Graph

        We can also use ``networkx``.

        .. doctest:: mod
            :options: +NORMALIZE_WHITESPACE

            >>> from gurobi_optimods.min_cost_flow import min_cost_flow_networkx
            >>> from gurobi_optimods import datasets
            >>> import networkx as nx
            >>>
            >>> edge_data, node_data = datasets.load_min_cost_flow2()
            >>> G = nx.from_pandas_edgelist(
            ...     edge_data.reset_index(), create_using=nx.DiGraph(), edge_attr=True
            ... )
            >>> for i, d in node_data.iterrows():
            ...     G.add_node(i, demand=d.demand, pos=(d.posx, d.posy))
            ...
            >>>
            >>> obj, sol = min_cost_flow_networkx(G) # doctest: +IGNORE_RESULT
            ...
            >>> edge_labels = {}
            >>> for i, j in G.edges():
            ...     f = 0.0
            ...     if (i, j) in sol:
            ...         f = sol[i, j]
            ...     edge_labels[i, j] = str((G.edges[i, j]["capacity"], G.edges[i, j]["cost"], f))
            ...
            >>> node_labels1 = {}
            >>> node_labels2 = {}
            >>> for n in G.nodes():
            ...     node_labels1[n] = f"{n}"
            ...     node_labels2[n] = f"{G.nodes[n]['demand']}"
            ...
            >>> color_map = ["red" if (i, j) in sol else "lightgrey" for (i, j) in G.edges()]
            >>> pos = nx.get_node_attributes(G, "pos")
            >>> nx.draw(G, with_labels=True, pos=pos, edge_color=color_map)
            >>> nx.draw_networkx_edge_labels(G, pos, edge_labels) # doctest: +IGNORE_RESULT
            >>> nx.draw_networkx_labels(G, pos, node_labels1) # doctest: +IGNORE_RESULT
            >>> for k, v in pos.items():
            ...     pos[k] = (v[0], v[1] + 0.07)
            ...
            >>> nx.draw_networkx_labels(G, pos, node_labels2, font_color="r") # doctest: +IGNORE_RESULT

The solution is shown in the following figure below. The edge labels denote the
edge capacity, cost and resulting flow: :math:`(B_{ij}, c_{ij}, x^*_{ij})`.
Edges with non-zero flow are highlighted in red. The demand for each vertex is
shown on top of the vertex in red. A negative demand indicates a supply vertex
whereas a positive demand indicates a consumer vertex.

.. image:: figures/min-cost-flow-result.png
  :width: 600
  :alt: Sample network.
