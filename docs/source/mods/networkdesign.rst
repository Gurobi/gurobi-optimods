.. This template should be copied to docs/source/mods/<mod_name>.rst

Network Design
==========

A little background on the proud history of mathprog in this field.

Also data science.

Problem Specification
---------------------

Give a brief overview of the problem being solved.


.. tabs::

    .. tab:: Graph Theory

        For a given graph :math:`G` with set of vertices :math:`V` and potential edges
        :math:`E` as well as a set of commodities :math:`K`.

        Each potential edge :math:`(i,j)\in E` has the following attributes:

        - flow cost: :math:`c_{ij}\in \mathbb{R}`;
        - fixed cost for building the arc: :math:`f_{ij} \in \mathbb{R}`;
        - and capacity: :math:`B_{ij}\in\mathbb{R}`.

        Each commodity :math:`k \in K` has the following attributes:

        - origin :math:`o_k \in V`;
        - destination: :math:`d_k \in V`;
        - and demand: :math:`D_k \in \mathbb{R}`


        Also, each vertex :math:`i\in V` has a demand :math:`d_i\in\mathbb{R}`.
        This value can be positive (requesting flow), negative (supplying
        flow), or 0.

        The problem can be stated as finding a the flow with minimal total cost
        such that:

        - the demand at each vertex is met;
        - and, the flow is capacity feasible.

    .. tab:: Optimization Model

        Let us define a set of continuous variables :math:`x_{ij}` to represent
        the amount of non-negative (:math:`\geq 0`) flow going through an edge
        :math:`(i,j)\in E`.


        The mathematical formulation can be stated as follows:

        .. math::


            \begin{alignat}{2}
            \min \quad &\sum_k \sum_{(i,j)\in A} W^k c_{ij} x^k_{ij} + \sum_{(i,j)\in A} f_{ij} y_{ij}\\
            \text{s.t.} \quad &\sum_{j\in \delta^+(i)} x^k_{ij} - \sum_{j\in \delta^-(i)} x^k_{ji} =\begin{cases}1 &\text{if }i=o_k\\-1 &\text{if }i=d_k\\0&\text{otherwise}\end{cases} \quad&\forall\ i\in N,\ k\in K\\
            &\sum_{k\in K} W^k x^k_{ij} \le C_{ij} y_{ij} & \forall\ (i,j)\in A\\
            & x^k_{ij}\geq 0,\ \ y_{ij}\in\{0,1\} &\forall\ (i,j)\in A,\ k\in K
            \end{alignat}

        Where :math:`\delta^+(\cdot)` (:math:`\delta^-(\cdot)`) denotes the
        outgoing (incoming) neighours.

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
Give examples of the various input data structures. These inputs should be fixed,
so use doctests where possible.

.. testsetup:: mod

    # Set pandas options for displaying dataframes, if needed
    import pandas as pd
    pd.options.display.max_rows = 10

.. tabs::

    .. tab:: ``availability``

        Give interpretation of input data.

        .. doctest:: mod
            :options: +NORMALIZE_WHITESPACE

            >>> from gurobi_optimods import datasets
            >>> data = datasets.load_workforce()
            >>> data.availability
               Worker      Shift
            0     Amy 2022-07-02
            1     Amy 2022-07-03
            2     Amy 2022-07-05
            3     Amy 2022-07-07
            4     Amy 2022-07-09
            ..    ...        ...
            67     Gu 2022-07-10
            68     Gu 2022-07-11
            69     Gu 2022-07-12
            70     Gu 2022-07-13
            71     Gu 2022-07-14
            <BLANKLINE>
            [72 rows x 2 columns]

        In the model, this corresponds to ...

    .. tab:: ``shift_requirements``

        Another bit of input data (perhaps a secondary table)

|

Code
----

Self contained code example to run the mod from an example dataset. Example
datasets should bd included in the ``gurobi_optimods.datasets`` module for
easy access by users.

.. testcode:: mod

    import pandas as pd

    from gurobi_optimods.datasets import load_mod_data
    from gurobi_optimods.mod import solve_mod


    data = load_mod_data()
    solution = solve_mod(data.table1, data.table2)

..  A snippet of the Gurobi log output here won't show in the rendered page,
    but serves as a doctest to make sure the code example runs. The ... lines
    are meaningful here, they will match anything in the output test.

.. testoutput:: mod
    :hide:

    ...
    Optimize a model with 14 rows, 72 columns, and 72 nonzeros
    ...
    Optimal objective
    ...

The model is solved as an LP/MIP/QP by Gurobi.

..  You can include the full Gurobi log output here for the curious reader.
    It will be visible as a collapsible section.

.. collapse:: View Gurobi Logs

    .. code-block:: text

        Gurobi Optimizer version 9.5.1 build v9.5.1rc2 (mac64[x86])
        Optimize a model with ...
        Best obj ... Best bound ...

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
