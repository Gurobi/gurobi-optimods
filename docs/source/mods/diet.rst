Stigler Diet Problem
====================

.. warning::
    This page is under construction. To put it bluntly, the diet mod has made
    some really poor diet choices and it should be ashamed of itself. But, as
    we know: garbage data in, garbage results out.

Thanks `Wikipedia <https://en.wikipedia.org/wiki/Stigler_diet>`_:

    For a moderately active man weighing 154 pounds, how much of each of 77 foods
    should be eaten on a daily basis so that the man’s intake of nine nutrients
    will be at least equal to the recommended dietary allowances (RDAs) suggested
    by the National Research Council in 1943, with the cost of the diet being minimal?

Problem Specification
---------------------

The classic diet problem involves a set of foods, each of which contain
varying amounts of each of several nutrients. Each food has a per-volume
cost. The problem requires us to minimize the total cost of a diet which
meets minimum requirements specified for each nutrient.

.. tabs::

    .. tab:: Description

        We have some data on various food items: their cost per unit,
        and amount per unit of various nutrients. We are also given
        a specification for a target diet: minimum and maximum quantities
        of each nutrient which may be consumed.

        Our goal is to find a diet (choice of food items and amounts of
        each) which satisfies the min and max limits for all required
        nutrients, for the lowest possible cost.

    .. tab:: Optimization Model

	    The model is defined over foods :math:`i \in F` and nutrients :math:`j
	    \in N`. Foods have cost per unit :math:`c_{i}` and nutrients per
	    unit :math:`n_{ij}`.

        .. math::

            \begin{alignat}{2}
            \min \quad        & \sum_{i \in F} c_{i} x_{i} \\
            \mbox{s.t.} \quad & l_{j} \le \sum_{i \in F} n_{ij} x_{i} \le u_{j} \,\, & \forall j \in N \\
            \end{alignat}

        Implementing this mathematical model using ``gurobipy-pandas`` is quite
        straightforward. The source for the mod can be found in the
        :ghsrc:`Github repository<src/gurobi_optimods/diet.py#L19>`.

The API uses pandas dataframes for each data table. Separate tables are needed
for foods and nutrients. We also need a linking table which specifies the nutrient
content of each food.

.. testsetup:: diet

    # Set pandas options for displaying dataframes, if needed
    import pandas as pd
    pd.options.display.max_rows = 10

.. tabs::

    .. tab:: ``categories``

        Data on dietary requirements is included as a dataframe.

        .. doctest:: diet
            :options: +NORMALIZE_WHITESPACE

            >>> from gurobi_optimods import datasets
            >>> data = datasets.load_diet()
            >>> data.categories
               category   min      max
            0  calories  1800   2200.0
            1   protein    91  10000.0
            2       fat     0     65.0
            3    sodium     0   1779.0

        The min and max columns correspond to :math:`l_j` and :math:`u_j` in the
        mathematical model.

    .. tab:: ``foods``

        Data on available foods is included as a dataframe.

        .. doctest:: diet
            :options: +NORMALIZE_WHITESPACE

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

	    The cost column corresponds to :math:`c_i` in the mathematical model.

    .. tab:: ``nutation``

        .. doctest:: diet
            :options: +NORMALIZE_WHITESPACE

            >>> data.nutrition_values
                category       food   value
            0   calories  hamburger   410.0
            1   calories    chicken   420.0
            2   calories    hot dog   560.0
            3   calories      fries   380.0
            4   calories   macaroni   320.0
            ..       ...        ...     ...
            31    sodium   macaroni   930.0
            32    sodium      pizza   820.0
            33    sodium      salad  1230.0
            34    sodium       milk   125.0
            35    sodium  ice cream   180.0
            <BLANKLINE>
            [36 rows x 3 columns]

        The value column corresponds to :math:`n_{ij}` in the mathematical model.

|

Code
----

Using the example data above, solve for the optimal diet.

.. testcode:: diet

    import pandas as pd

    from gurobi_optimods.datasets import load_diet
    from gurobi_optimods.diet import solve_diet_problem


    data = load_diet()
    solution = solve_diet_problem(
        categories=data.categories,
        foods=data.foods,
        values=data.nutrition_values,
    )

.. testoutput:: diet
    :hide:

    ...
    Optimize a model with 8 rows, 9 columns and 72 nonzeros
    ...

Gurobi solves this now-simple linear programming model with ease. But think back
to the days of Dantzig where hundreds of engineers performed simplex iterations
using slide rules and sextants.

.. collapse:: View Gurobi Logs

    .. code-block:: text

        Gurobi Optimizer version 10.0.0 build v10.0.0rc2 (mac64[x86])

        CPU model: Intel(R) Core(TM) i5-1038NG7 CPU @ 2.00GHz
        Thread count: 4 physical cores, 8 logical processors, using up to 8 threads

        Optimize a model with 8 rows, 9 columns and 72 nonzeros
        Model fingerprint: 0x4ec4fbc2
        Coefficient statistics:
        Matrix range     [2e+00, 2e+03]
        Objective range  [9e-01, 3e+00]
        Bounds range     [0e+00, 0e+00]
        RHS range        [6e+01, 1e+04]
        Presolve removed 4 rows and 0 columns
        Presolve time: 0.00s
        Presolved: 4 rows, 10 columns, 37 nonzeros

        Iteration    Objective       Primal Inf.    Dual Inf.      Time
            0    0.0000000e+00   1.472500e+02   0.000000e+00      0s
            4    1.1828861e+01   0.000000e+00   0.000000e+00      0s

        Solved in 4 iterations and 0.00 seconds (0.00 work units)
        Optimal objective  1.182886111e+01

|

Solution
--------

Here's the result! We display as a simple pandas series for now. The result
object also contains the total cost of the menu as an attribute.

.. doctest:: diet
    :options: +NORMALIZE_WHITESPACE

    >>> solution.menu.round(2)
    food
    hamburger    0.60
    chicken      0.00
    hot dog      0.00
    fries        0.00
    macaroni     0.00
    pizza        0.00
    salad        0.00
    milk         6.97
    ice cream    2.59
    Name: quantity, dtype: float64
    >>> round(solution.total_cost, 2)
    11.83
