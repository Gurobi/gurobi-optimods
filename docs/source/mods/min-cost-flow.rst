Minimum-Cost Flow
=================

Minimum-cost flow problems are defined on a graph where the goal is to route
a certain amount of flow in the cheapest way. It is a fundamental flow problem
as many other graph problems can be modelled using this framework, for example,
the shortest-path, maximum flow, or matching problems.

The first algorithm to solve this problem was proposed by Dantzig
:footcite:p:`dantzig1951application`, :footcite:p:`dantzig1963linear`, called
the `network simplex` (NS) algorithm. Other methods have been proposed since
then but NS remains one the most efficient approaches. Competitive methods on
larger networks include cost-scaling methods (e.g. variants of the push-relabel
algorithm by :footcite:t:`goldberg1990finding`). For a more detailed comparison
see, for example, :footcite:t:`kovacs2015minimum`.

Problem Specification
---------------------

For a given graph :math:`G` with set of vertices :math:`V` and edges
:math:`E`. Each edge :math:`(i,j)\in E` has the following attributes:

- cost: :math:`c_{ij}\in \mathbb{R}`;
- and capacity: :math:`B_{ij}\in\mathbb{R}`.

Each vertex :math:`i\in V` has a demand :math:`d_i\in\mathbb{R}` which can be
positive (requesting flow), negative (supplying flow), or zero (a transshipment
node).

The problem can be stated as finding the flow with minimal total cost
such that:

- the demand at each vertex is met exactly; and
- the flow is capacity feasible.

.. dropdown:: Background: Optimization Model

    This Mod is implemented by formulating a Linear Program (LP) and solving it
    using Gurobi. Let us define a set of continuous variables :math:`x_{ij}` to
    represent the amount of non-negative (:math:`\geq 0`) flow going through an
    edge :math:`(i,j)\in E`.

    The formulation can be stated as follows:

    .. math::

        \begin{alignat}{2}
            \min \quad        & \sum_{(i, j) \in E} c_{ij} x_{ij} \\
            \mbox{s.t.} \quad & \sum_{j \in \delta^+(i)} x_{ij} - \sum_{j \in \delta^-(i)} x_{ji} = d_i & \forall i \in V \\
                            & 0 \leq x_{ij} \le B_{ij} & \forall (i, j) \in E \\
        \end{alignat}

    Where :math:`\delta^+(\cdot)` (:math:`\delta^-(\cdot)`) denotes the
    outgoing (incoming) neighbours.

    The objective minimises the total cost over all edges.

    The first constraints ensure flow balance for all vertices. That is, for
    a given node, the incoming flow (sum over all incoming edges to this
    node) minus the outgoing flow (sum over all outgoing edges from this
    node) is equal to the demand. Clearly, in the case when the demand is 0,
    the outgoing flow must be equal to the incoming flow. When the demand is
    negative, this node can supply flow to the network (outgoing term is
    larger), and conversely when the demand is negative, this node can
    request flow from the network (incoming term is larger).

    The last constraints ensure non-negativity of the variables and that the
    capacity per edge is not exceeded.


Code and Inputs
---------------

For this mod, one can use input graphs of different types:

* pandas: using a ``pd.DataFrame``;
* NetworkX: using a ``nx.DiGraph`` or ``nx.Graph``;
* SciPy.sparse: using some ``sp.sparray`` matrices and NumPy's ``np.ndarray``.

An example of these inputs with their respective requirements is shown below.

.. tabs::

  .. group-tab:: scipy.sparse

      .. doctest:: load_graph_scipy
          :options: +NORMALIZE_WHITESPACE

          >>> from gurobi_optimods import datasets
          >>> G, capacities, cost, demands = datasets.simple_graph_scipy()
          >>> G
          <5x6 sparse matrix of type '<class 'numpy.int64'>'
                  with 7 stored elements in COOrdinate format>
          >>> print(G)
            (0, 1)        1
            (0, 2)        1
            (1, 3)        1
            (2, 3)        1
            (2, 4)        1
            (3, 5)        1
            (4, 5)        1
          >>> print(capacities)
            (0, 1)        2
            (0, 2)        2
            (1, 3)        1
            (2, 3)        1
            (2, 4)        2
            (3, 5)        2
            (4, 5)        2
          >>> print(cost)
            (0, 1)        9
            (0, 2)        7
            (1, 3)        1
            (2, 3)        10
            (2, 4)        6
            (3, 5)        1
            (4, 5)        1
          >>> print(demands)
          [-2  0 -1  1  0  2]

      Three separate sparse matrices including the adjacency matrix, edge
      capacity and cost, and a single array with the demands per node.

  .. group-tab:: networkx

      .. doctest:: load_graph_networkx
          :options: +NORMALIZE_WHITESPACE

          >>> from gurobi_optimods import datasets
          >>> G = datasets.simple_graph_networkx()
          >>> for e in G.edges(data=True):
          ...     print(e)
          ...
          (0, 1, {'capacity': 2, 'cost': 9})
          (0, 2, {'capacity': 2, 'cost': 7})
          (1, 3, {'capacity': 1, 'cost': 1})
          (2, 3, {'capacity': 1, 'cost': 10})
          (2, 4, {'capacity': 2, 'cost': 6})
          (3, 5, {'capacity': 2, 'cost': 1})
          (4, 5, {'capacity': 2, 'cost': 1})
          >>> for n in G.nodes(data=True):
          ...     print(n)
          ...
          (0, {'demand': -2})
          (1, {'demand': 0})
          (2, {'demand': -1})
          (3, {'demand': 1})
          (4, {'demand': 0})
          (5, {'demand': 2})

      Edges have attributes ``capacity`` and ``cost`` and nodes have
      attributes ``demand``.

  .. group-tab:: pandas

      .. doctest:: load_graph
          :options: +NORMALIZE_WHITESPACE

          >>> from gurobi_optimods import datasets
          >>> edge_data, node_data = datasets.simple_graph_pandas()
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
             demand
          0      -2
          1       0
          2      -1
          3       1
          4       0
          5       2

      The ``edge_data`` DataFrame is indexed by ``source`` and ``target``
      nodes and contains columns labelled ``capacity`` and ``cost`` with the
      edge attributes.

      The ``node_data`` DataFrame is indexed by node and contains columns
      labelled ``demand``.

      We assume that nodes labels are integers from :math:`0,\dots,|V|-1`.

Solution
--------

Depending on the input of choice, the solution also comes with different
formats.

.. tabs::

  .. group-tab:: scipy.sparse

      .. doctest:: min_cost_flow_networkx
          :options: +NORMALIZE_WHITESPACE

          >>> from gurobi_optimods import datasets
          >>> from gurobi_optimods.min_cost_flow import min_cost_flow_scipy
          >>> G, capacities, cost, demands = datasets.simple_graph_scipy()
          >>> obj, sol = min_cost_flow_scipy(G, capacities, cost, demands, verbose=False)
          >>> obj
          31.0
          >>> sol
          <5x6 sparse matrix of type '<class 'numpy.float64'>'
                  with 5 stored elements in COOrdinate format>
          >>> print(sol)
            (0, 1)        1.0
            (0, 2)        1.0
            (1, 3)        1.0
            (2, 4)        2.0
            (4, 5)        2.0

      The ``min_cost_flow_scipy`` function returns the cost of the solution as
      well as a ``sp.sparray`` with the edges where the data is the amount of
      non-zero flow in the solution.

  .. group-tab:: networkx

      .. doctest:: min_cost_flow_networkx
          :options: +NORMALIZE_WHITESPACE

          >>> from gurobi_optimods import datasets
          >>> from gurobi_optimods.min_cost_flow import min_cost_flow_networkx
          >>> G = datasets.simple_graph_networkx()
          >>> obj, sol = min_cost_flow_networkx(G, verbose=False)
          >>> obj
          31.0
          >>> list(sol.edges(data=True))
          [(0, 1, {'flow': 1.0}), (0, 2, {'flow': 1.0}), (1, 3, {'flow': 1.0}), (2, 4, {'flow': 2.0}), (4, 5, {'flow': 2.0})]

      The ``min_cost_flow_networkx`` function returns the cost of the solution
      as well as a dictionary indexed by edge with the non-zero flow.

  .. group-tab:: pandas

      .. doctest:: min_cost_flow
          :options: +NORMALIZE_WHITESPACE

          >>> from gurobi_optimods import datasets
          >>> from gurobi_optimods.min_cost_flow import min_cost_flow_pandas
          >>> edge_data, node_data = datasets.simple_graph_pandas()
          >>> obj, sol = min_cost_flow_pandas(edge_data, node_data, verbose=False)
          >>> obj
          31.0
          >>> sol
          source  target
          0       1         1.0
                  2         1.0
          1       3         1.0
          2       3         0.0
                  4         2.0
          3       5         0.0
          4       5         2.0
          Name: flow, dtype: float64

      The ``min_cost_flow_pandas`` function returns the cost of the solution as well
      as ``pd.Series`` with the flow per edge. Similarly as the input
      DataFrame the resulting series is indexed by ``source`` and ``target``.

The solution for this example is shown in the figure below. The edge labels
denote the edge capacity, cost and resulting flow: :math:`(B_{ij}, c_{ij},
x^*_{ij})`. Edges with non-zero flow are highlighted in red. Also the demand for
each vertex is shown on top of the vertex in red.

.. image:: figures/min-cost-flow-result.png
  :width: 600
  :alt: Sample network.

In all these cases, the model is solved as an LP by Gurobi (typically using the
NS algorithm).

.. footbibliography::
