Workforce Scheduling
====================

A little background on the proud history of mathprog in workforce scheduling.

Also data science.

Todo:

- Add a ceiling on the number of shifts fper worker (worker data or fixed value)
- Add an optional rolling time window constraint (rotating roster)

Note that simplified versions of this model can be solved as graph matching
problems:

* If we ignore all worker preferences and impose no limits on the number of
  shifts assigned to any given worker, then the problem of assigning available
  workers to shifts can be solved by maximum cardinality bipartite matching.
* If we include preferences, the problem can be solved using weighted

Once side constraints (e.g. total shifts per worker) get involved, we need to
move away from these simplified models. However the same basic formulation
approach of declaring a binary variable for every possible worker-shift
assignment can be used throughout.

Problem Specification
---------------------

Consider a service business, like a restaurant, that develops its workforce plans for the next two weeks (considering a 7-day week). The service requires only one set of skills. There are a number of employed workers with the same set of skills and with identical productivity that are available to work on some of the days during the two-week planning horizon. There is only one shift per workday. Each shift may have different resource (worker) requirements on each workday. The service business wants to minimize the cost of covering all shifts with available workers, based on their rates of pay and availability.

.. tabs::

    .. tab:: Data Specification

        The workforce scheduling model takes the following inputs:

        * The ``availability`` dataframe has three columns: ``Worker``, ``Shift``,
          and ``Preference``. Each row in dataframe specifies that the given worker
          is available to work the given shift.
        * The ``shift_requirements`` dataframe has two columns: ``Shift`` and
          ``Required``. Each row specifies the number of workers required for a
          given shift. There should be one row for every unique worker in
          ``availability["Workers"]``.

        When ``solve_workforce_scheduling`` is called, a model is formulated and
        solved immediately using Gurobi. Workers will be assigned only to shifts
        they are available for, in such a way that all requirements are covered while
        total sum of worker preference scores is maximised.

        The returned assignment dataframe is a subset of the availability dataframe,
        with the same columns. Each row specifies that the given worker has been
        assigned to the given shift.

    .. tab:: Optimization Model

        Set of :math:`S` shifts to cover using set of workers :math:`W`. Workers
        :math:`w \in W_{s} \subseteq W` are available to work a given shift `s`,
        and have a preference :math:`p_{ws}` for each assigned shift. Shift :math:`s`
        requires :math:`r_{s}` workers assigned. The model is defined on variables
        :math:`x_{ws}` such that

        .. math::

            x_{ws} = \begin{cases}
                1 & \text{if worker w is given shift s} \\
                0 & \text{otherwise.} \\
            \end{cases}

        The mathematical model is then expressed as:

        .. math::

            \begin{alignat}{2}
            \max \quad        & \sum_{s \in S} \sum_{w \in W_{s}} p_{ws} x_{ws} \\
            \mbox{s.t.} \quad & \sum_{w \in W_{s}} x_{ws} = r_{s} & \forall s \in S \\
                              & 0 \le x_{ws} \le 1 & \forall s \in S, w \in W_{s} \\
            \end{alignat}

        The objective computes the total cost of all shift assignments based on
        worker pay rates. The constraint ensures all shifts are assigned the
        required number of workers.

        Technically, :math:`x_{ws}` should be binary, but we're in magical
        assignment land here. Read more about why assignment problems and network
        flows are the bees knees `here <www.gurobi.com>`_. Note that if the model is
        modified with additional constraints, this elegant property may no longer hold.

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

    .. tab:: ``preferences``

        Amy is available for a shift on Tuesday 2nd, etc, etc

        .. doctest:: workforce
            :options: +NORMALIZE_WHITESPACE

            >>> from gurobi_optimods import datasets
            >>> data = datasets.load_workforce()
            >>> data.preferences
               Worker      Shift  Preference
            0     Amy 2022-07-02         3.0
            1     Amy 2022-07-03         3.0
            2     Amy 2022-07-05         3.0
            3     Amy 2022-07-07         3.0
            4     Amy 2022-07-09         3.0
            ..    ...        ...         ...
            67     Gu 2022-07-10         2.0
            68     Gu 2022-07-11         2.0
            69     Gu 2022-07-12         2.0
            70     Gu 2022-07-13         2.0
            71     Gu 2022-07-14         2.0
            <BLANKLINE>
            [72 rows x 3 columns]

        In the mathematical model, this models the set :math:`\lbrace (w, s) \mid s \in S, w \in W_s \rbrace` and preference values :math:`p_{ws}`.

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
        availability=data.preferences,
        shift_requirements=data.shift_requirements,
    )

.. testoutput:: workforce
    :hide:

    ...
    Optimize a model with 14 rows, 72 columns and 72 nonzeros
    ...
    Best objective 1.960000000000e+02, best bound 1.960000000000e+02, gap 0.0000%

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
       Worker      Shift  Preference
    0     Amy 2022-07-03         3.0
    1     Amy 2022-07-05         3.0
    2     Amy 2022-07-07         3.0
    3     Amy 2022-07-10         3.0
    4     Amy 2022-07-11         3.0
    ..    ...        ...         ...
    47     Gu 2022-07-05         2.0
    48     Gu 2022-07-06         2.0
    49     Gu 2022-07-07         2.0
    50     Gu 2022-07-12         2.0
    51     Gu 2022-07-13         2.0
    <BLANKLINE>
    [52 rows x 3 columns]

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

Further Constraints
-------------------

In order to satisfy employee agreements, we implement some simple constraints
for this two-week period. Namely, each worker will have a minimum number of
shifts they are entitled to, and will not be allocated more shifts than a given
maximum.

Some additional data needs to be provided to achieve this, in the form of a
``worker_data`` dataframe ...

Rolling Rosters
---------------

Let's now assume that the above requirements must be satisfied on a rolling
basis. So, instead of providing upper and lower limits on the number of shifts
assigned in this 14 day period, we enforce similar limits on a rolling basis.

This is a stricter constraint: in `any` 14 day window defined by starting at
a given day in the roster, enforce the upper and lower limits.

The constraint to add is no more complex, we just need to implement the code.
