.. This template should be copied to docs/source/mods/<mod_name>.rst

Stigler Diet Problem
====================

Thanks `Wikipedia <https://en.wikipedia.org/wiki/Stigler_diet>`_!

For a moderately active man weighing 154 pounds, how much of each of 77 foods
should be eaten on a daily basis so that the man’s intake of nine nutrients
will be at least equal to the recommended dietary allowances (RDAs) suggested
by the National Research Council in 1943, with the cost of the diet being minimal?

Problem Specification
---------------------

Give a brief overview of the problem being solved.

.. tabs::

    .. tab:: Domain-Specific Description

        We have some data on various food items: their cost per unit,
        and amount per unit of various nutrients. We are also given
        a specification for a target diet: minimum and maximum quantities
        of each nutrient which may be consumed.

        Our goal is to find a diet (choice of food items and amounts of
        each) which satisfies the min and max limits for all required
        nutrients, for the lowest possible cost.

    .. tab:: Optimization Model

	The model is defined over foods :math:`i \in F` and nutrients :math:`j \in N`.
        Foods have cost per unit :math:`c_{i}` and nutrients per unit :math:`n_{ij}`.

        .. math::

            \begin{alignat}{2}
            \min \quad        & \sum_{i \in F} c_{i} x_{i} \\
            \mbox{s.t.} \quad & l_{j} \le \sum_{i \in F} n_{ij} x_{i} \le u_{j} \,\, & \forall j \in N \\
            \end{alignat}

Give examples of the various input data structures. These inputs should be fixed,
so use doctests where possible.

.. testsetup:: mod

    # Set pandas options for displaying dataframes, if needed
    import pandas as pd
    pd.options.display.max_rows = 10

.. tabs::

    .. tab:: ``categories``

        Give interpretation of input data.

        .. doctest:: mod
            :options: +NORMALIZE_WHITESPACE

            >>> from gurobi_optimods import datasets
            >>> data = datasets.load_diet()
            >>> data.categories
               category   min      max
            0  calories  1800   2200.0
            1   protein    91  10000.0
            2       fat     0     65.0
            3    sodium     0   1779.0

        The min and max columns correspond to :math:`l_j` and :math:`u_j`.

    .. tab:: ``foods``

        Another bit of input data (perhaps a secondary table)

            >>> from gurobi_optimods import datasets
            >>> data = datasets.load_diet()
            >>> data.foods
                    food  cost
            0  hamburger  2.49
            1    chicken  2.89
            2    hot dog  1.50
            3      fries  1.89
            4   macaroni  2.09
            5      pizza  1.99
            6      salad  2.49
            7       milk  0.89
            8  ice cream  1.59

	The cost column corresponds to :math:`c_i`.

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
