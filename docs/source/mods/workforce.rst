Workforce Scheduling
====================

A little background on the proud history of mathprog in workforce scheduling.

Also data science.

- Convert cost minimization to preference maximisation
- Add an optional rolling time window constraint (rotating roster?)
- Mention bipartite and weighted matching as simple cases of the problem. Here we add more interesting practical constraints

Problem Specification
---------------------

Consider a service business, like a restaurant, that develops its workforce plans for the next two weeks (considering a 7-day week). The service requires only one set of skills. There are a number of employed workers with the same set of skills and with identical productivity that are available to work on some of the days during the two-week planning horizon. There is only one shift per workday. Each shift may have different resource (worker) requirements on each workday. The service business wants to minimize the cost of covering all shifts with available workers, based on their rates of pay and availability.

.. tabs::

    .. tab:: Data Specification

        The workforce scheduling model takes the following inputs:

        * The ``availability`` dataframe has two columns: ``Worker`` and ``Shift``. Each row in dataframe specifies that the given worker is available to work the given shift.
        * The ``shift_requirements`` dataframe has two columns: ``Shift`` and ``Required``. Each row specifies the number of workers required for a given shift. There should be one row for every unique worker in ``availability["Workers"]``.
        * The ``pay_rates`` dataframe has two columns: ``Worker`` and ``PayRate``. Each row specifies the pay per shift worked for that worker. There should be one row for every unique shift in ``availability["Shift"]``.

        When ``solve_workforce_scheduling`` is called, a model is formulated and solved immediately using Gurobi. Workers will be assigned only to shifts they are available for, in such a way that all requirements are covered while total cost of covering all shifts is minimised.

        The returned assignment dataframe is a subset of the availability dataframe, with the same columns. Each row specifies that the given worker has been assigned the given shift.

    .. tab:: Optimization Model

        Set of :math:`S` shifts to cover using set of workers :math:`W`. Workers :math:`w \in W_{s} \subseteq W` are available to work a given shift `s`, and are paid an amount :math:`c_{w}` for each assigned shift. Shift :math:`s` requires :math:`r_{s}` workers assigned. The model is defined on variables :math:`x_{ws}` such that

        .. math::

            x_{ws} = \begin{cases}
                1 & \text{if worker w is given shift s} \\
                0 & \text{otherwise.} \\
            \end{cases}

        The mathematical model is then expressed as:

        .. math::

            \begin{alignat}{2}
            \min \quad        & \sum_{s \in S} \sum_{w \in W_{s}} c_{w} x_{ws} \\
            \mbox{s.t.} \quad & \sum_{w \in W_{s}} x_{ws} = r_{s} & \forall s \in S \\
                              & 0 \le x_{ws} \le 1 & \forall s \in S, w \in W_{s} \\
            \end{alignat}

        The objective computes the total cost of all shift assignments based on worker pay rates.
        The constraint ensures all shifts are assigned the required number of workers.

        Technically, :math:`x_{ws}` should be binary, but we're in magical assignment land here.
        Read more about why assignment problems and network flows are the bees knees `here <www.gurobi.com>`_.
        Note that if the model is modified with additional constraints, this elegant property may no longer hold.

Data examples

.. a few re-use hacks for doctests

.. testsetup:: workforce

    # Set pandas options
    import pandas as pd
    pd.options.display.max_rows = 10

.. Maybe the example paths should be found in a datasets module
.. similar to sklearn. We could included proccessing code to
.. read from csv and avoid the feather dependency that way.

.. tabs::

    .. tab:: ``availability``

        Amy is available for a shift on Tuesday 2nd, etc, etc

        .. doctest:: workforce
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

        In the mathematical model, this models the set :math:`\lbrace (w, s) \mid s \in S, w \in W_s \rbrace`.

    .. tab:: ``shift_requirements``

        Shift on Monday 1st requires 3 workers, etc, etc

        .. doctest:: workforce
            :options: +NORMALIZE_WHITESPACE

            >>> from gurobi_optimods import datasets
            >>> data = datasets.load_workforce()
            >>> data.shift_requirements
                    Shift  Required
            0  2022-07-01         3
            1  2022-07-02         2
            2  2022-07-03         4
            3  2022-07-04         2
            4  2022-07-05         5
            ..        ...       ...
            9  2022-07-10         3
            10 2022-07-11         4
            11 2022-07-12         5
            12 2022-07-13         7
            13 2022-07-14         5
            <BLANKLINE>
            [14 rows x 2 columns]

        In the mathematical model, this models the values :math:`r_s`.

    .. tab:: ``pay_rates``

        Bob is the most expensive worker ...

        .. doctest:: workforce
            :options: +NORMALIZE_WHITESPACE

            >>> from gurobi_optimods import datasets
            >>> data = datasets.load_workforce()
            >>> data.pay_rates
              Worker  PayRate
            0    Amy       10
            1    Bob       12
            2  Cathy       10
            3    Dan        8
            4     Ed        8
            5   Fred        9
            6     Gu       11

        In the mathematical model, this models the values :math:`c_w`.

|

Code
----

Show the code required to run the mod. Users interact with the 'solver' by passing dataframes to a given spec and receiving a dataframe as output.

.. testcode:: workforce

    import pandas as pd
    pd.options.display.max_rows = 15

    from gurobi_optimods.datasets import load_workforce
    from gurobi_optimods.workforce import solve_workforce_scheduling


    # Load example data.
    data = load_workforce()

    # Get winning results.
    assigned_shifts = solve_workforce_scheduling(
        availability=data.availability,
        shift_requirements=data.shift_requirements,
        pay_rates=data.pay_rates,
    )

.. testoutput:: workforce
    :hide:

    ...
    Optimize a model with 14 rows, 72 columns and 72 nonzeros
    ...
    Optimal objective  4.800000000e+02

The model is solved as a linear program by Gurobi.

.. collapse:: View Gurobi logs

    .. code-block:: text

        Gurobi Optimizer version 9.5.1 build v9.5.1rc2 (mac64[x86])
        Thread count: 4 physical cores, 8 logical processors, using up to 8 threads
        Optimize a model with 14 rows, 72 columns and 72 nonzeros
        Model fingerprint: 0x494be3a7
        Coefficient statistics:
        Matrix range     [1e+00, 1e+00]
        Objective range  [8e+00, 1e+01]
        Bounds range     [1e+00, 1e+00]
        RHS range        [2e+00, 7e+00]
        Presolve removed 14 rows and 72 columns
        Presolve time: 0.00s
        Presolve: All rows and columns removed
        Iteration    Objective       Primal Inf.    Dual Inf.      Time
            0    4.8000000e+02   0.000000e+00   1.480000e+02      0s
        Extra simplex iterations after uncrush: 5
            5    4.8000000e+02   0.000000e+00   0.000000e+00      0s

        Solved in 5 iterations and 0.00 seconds (0.00 work units)
        Optimal objective  4.800000000e+02

|

Solution
--------

Solution is a selection of shift assignments. The returned dataframe is just a
subset of the availability dataframe, so we can transform the results using
normal pandas code (no gurobipy interaction).

.. doctest:: workforce
    :options: +NORMALIZE_WHITESPACE

    >>> assigned_shifts
    Worker      Shift
    0     Amy 2022-07-03
    1     Amy 2022-07-05
    2     Amy 2022-07-07
    3     Amy 2022-07-10
    4     Amy 2022-07-11
    ..    ...        ...
    47     Gu 2022-07-05
    48     Gu 2022-07-06
    49     Gu 2022-07-07
    50     Gu 2022-07-12
    51     Gu 2022-07-13
    <BLANKLINE>
    [52 rows x 2 columns]

Further transform

.. doctest:: workforce
    :options: +NORMALIZE_WHITESPACE

    >>> shifts_table = pd.pivot_table(
    ...     assigned_shifts.assign(value=1),
    ...     values="value",
    ...     index="Shift",
    ...     columns="Worker",
    ...     fill_value="-",
    ... ).replace({1.0: "Y"})
    >>> shifts_table
    Worker     Amy Bob Cathy Dan Ed Fred Gu
    Shift
    2022-07-01   -   -     -   -  Y    Y  Y
    2022-07-02   -   -     -   Y  Y    -  -
    2022-07-03   Y   -     -   Y  Y    Y  -
    2022-07-04   -   -     Y   -  Y    -  -
    2022-07-05   Y   -     Y   Y  Y    -  Y
    2022-07-06   -   Y     -   Y  -    Y  Y
    2022-07-07   Y   -     Y   -  Y    -  Y
    2022-07-08   -   -     -   Y  Y    -  -
    2022-07-09   -   -     -   Y  Y    -  -
    2022-07-10   Y   -     Y   Y  -    -  -
    2022-07-11   Y   -     Y   Y  Y    -  -
    2022-07-12   Y   -     Y   Y  -    Y  Y
    2022-07-13   Y   Y     Y   Y  Y    Y  Y
    2022-07-14   Y   -     Y   Y  Y    Y  -
