Workforce Scheduling
====================

Workforce scheduling is an extremely widely-used application of optimization in
practice. It involves balancing many competing concerns, including worker
availability and preferences, shift coverage requirements, conditions on
consecutive shifts or rest breaks, and so on. Implementation can become quite
involved as worker requirements and entitlements become more complex.

This mod implements a relatively simple case. Workers provide their availability
and preferences, while rest requirements and work entitlements are handled
through lower and upper limits on the number of shifts a worker is rostered on
for in a given period. The scheduler aims to maximize satisfaction by finding a
feasible roster which maximizes the sum of preference scores.

Problem Specification
---------------------

Consider a service business, such as a restaurant, that develops its roster for
a two week period. The service requires only one set of skills. There are a
number of employed workers with the same set of skills and with identical
productivity that are available to work on some of the days during the two-week
planning horizon. There is only one shift per workday, however each shift may
require a different number of workers. The business requests preferences from
all workers for shifts they are available for, and aims to maximize the sum of
preference scores of assigned shifts as a proxy for worker happiness.

.. tabs::

    .. tab:: Data Specification

        The workforce scheduling model takes the following three dataframes as
        input:

        * The ``preferences`` dataframe has three columns: ``Worker``, ``Shift``,
          and ``Preference``. Each row in dataframe specifies that the given worker
          is available to work the given shift.
        * The ``shift_requirements`` dataframe has two columns: ``Shift`` and
          ``Required``. Each row specifies the number of workers required for a
          given shift. There must be one row for every unique shift in
          ``preferences["Shift"]``.
        * The ``worker_limits`` dataframe has three columns: ``Worker``,
          ``MinShifts``, and ``MaxShifts``. Each row specifies the minimum and
          maximum number of shifts the given worker may be assigned in the
          schedule. There must be one row for every unique worker in
          ``preferences["Worker"]``.

        When ``solve_workforce_scheduling`` is called, a model is formulated and
        solved immediately using Gurobi. Workers will be assigned only to shifts
        they are available for, in such a way that all requirements are covered,
        minimum and maximum shift numbers are respected, and the total sum of
        worker preference scores is maximised.

        The returned assignment dataframe is a subset of the preferences
        dataframe, with the same columns. A row in the returned dataframe
        specifies that the given worker has been assigned to the given shift.

    .. tab:: Optimization Model

        The set of shifts :math:`S` is to be covered using the set of workers
        :math:`W`. Workers :math:`w \in W_{s} \subseteq W` are available to work
        a given shift `s`, and have a preference :math:`p_{ws}` for each
        assigned shift. Shift :math:`s` requires :math:`r_{s}` workers assigned.
        Each worker must be assigned between :math:`l_{w}` and :math:`u_{w}`
        shifts in total. The model is defined on binary variables :math:`x_{ws}`
        which satisfy the condition

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
                              & l_{w} \le \sum_{s \in S} x_{ws} \le u_{w} & \forall w \in W \\
                              & x_{ws} \in \lbrace 0, 1 \rbrace & \forall s \in S, w \in W_{s} \\
            \end{alignat}

        The objective computes the total cost of all shift assignments based on
        worker pay rates. The first constraint ensures that all shifts are
        assigned the required number of workers, while the second constraint
        ensures workers are assigned to an acceptable number of shifts.

All input data is given as pandas dataframes, following the layout described in
the Data Specification above. The tabs below show example data for each frame.

.. testsetup:: workforce

    # Set pandas options
    import pandas as pd
    pd.options.display.max_rows = 10

.. tabs::

    .. tab:: ``preferences``

        The following example table lists worker availability and preferences.
        For example, Amy is available on July 2nd, 3rd, 5th, and so on, with a
        stronger preference to be assigned the shift on the 5th.

        .. doctest:: workforce
            :options: +NORMALIZE_WHITESPACE

            >>> from gurobi_optimods import datasets
            >>> data = datasets.load_workforce()
            >>> data.preferences
               Worker      Shift  Preference
            0     Amy 2022-07-02         2.0
            1     Amy 2022-07-03         2.0
            2     Amy 2022-07-05         5.0
            3     Amy 2022-07-07         3.0
            4     Amy 2022-07-09         2.0
            ..    ...        ...         ...
            67     Gu 2022-07-10         4.0
            68     Gu 2022-07-11         5.0
            69     Gu 2022-07-12         2.0
            70     Gu 2022-07-13         4.0
            71     Gu 2022-07-14         3.0
            <BLANKLINE>
            [72 rows x 3 columns]

        In the mathematical model, the worker-shift pairings model the set
        :math:`\lbrace (w, s) \mid s \in S, w \in W_s \rbrace` and the
        preference column provides values :math:`p_{ws}`.

    .. tab:: ``shift_requirements``

        The following example table lists the number of workers required for
        each shift.

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

        In the mathematical model, this table provides the values :math:`r_s`.

    .. tab:: ``worker_limits``

        The following example table lists the minimum and maximum number of
        shifts in the planning period which each worker is entitled to.

        .. doctest:: workforce
            :options: +NORMALIZE_WHITESPACE

            >>> from gurobi_optimods import datasets
            >>> data = datasets.load_workforce()
            >>> data.worker_limits
              Worker  MinShifts  MaxShifts
            0    Amy          6          8
            1    Bob          6          7
            2  Cathy          6          8
            3    Dan          5          8
            4     Ed          6          8
            5   Fred          5          8
            6     Gu          6          8

        In the mathematical model, this table provides the values :math:`l_w`
        and :math:`u_w`.

Solving a Model
---------------

The example code below solves the workforce scheduling problem for a simple
example dataset comprising seven workers covering daily shifts over a two week
period.

.. testcode:: workforce

    from gurobi_optimods.datasets import load_workforce
    from gurobi_optimods.workforce import solve_workforce_scheduling

    # Load example data
    data = load_workforce()

    # Solve the mod, get back a schedule
    assigned_shifts = solve_workforce_scheduling(
        preferences=data.preferences,
        shift_requirements=data.shift_requirements,
        worker_limits=data.worker_limits,
    )

.. testoutput:: workforce
    :hide:

    ...
    Optimize a model with 28 rows, 72 columns and 216 nonzeros
    ...
    Best objective 1.850000000000e+02, best bound 1.850000000000e+02, gap 0.0000%

.. collapse:: View Gurobi logs for solving this example

    .. code-block:: text

        Gurobi Optimizer version 10.0.1 build v10.0.1rc0 (mac64[x86])

        CPU model: Intel(R) Core(TM) i5-1038NG7 CPU @ 2.00GHz
        Thread count: 4 physical cores, 8 logical processors, using up to 8 threads

        Optimize a model with 28 rows, 72 columns and 216 nonzeros
        Model fingerprint: 0x595b329f
        Variable types: 0 continuous, 72 integer (72 binary)
        Coefficient statistics:
          Matrix range     [1e+00, 1e+00]
          Objective range  [1e+00, 5e+00]
          Bounds range     [1e+00, 1e+00]
          RHS range        [2e+00, 8e+00]
        Found heuristic solution: objective 170.0000000
        Presolve removed 6 rows and 22 columns
        Presolve time: 0.00s
        Presolved: 22 rows, 50 columns, 145 nonzeros
        Variable types: 0 continuous, 50 integer (50 binary)
        Found heuristic solution: objective 177.0000000

        Root relaxation: objective 1.850000e+02, 24 iterations, 0.00 seconds (0.00 work units)

            Nodes    |    Current Node    |     Objective Bounds      |     Work
         Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

        *    0     0               0     185.0000000  185.00000  0.00%     -    0s

        Explored 1 nodes (24 simplex iterations) in 0.00 seconds (0.00 work units)
        Thread count was 8 (of 8 available processors)

        Solution count 3: 185 177 170

        Optimal solution found (tolerance 1.00e-04)
        Best objective 1.850000000000e+02, best bound 1.850000000000e+02, gap 0.0000%

|

Inspecting the Solution
-----------------------

The solution to this workforce scheduling problem is a selection of shift
assignments. The returned dataframe is a subset of the original preferences
dataframe.

.. doctest:: workforce
    :options: +NORMALIZE_WHITESPACE

    >>> assigned_shifts
       Worker      Shift  Preference
    0     Amy 2022-07-03         2.0
    1     Amy 2022-07-05         5.0
    2     Amy 2022-07-07         3.0
    3     Amy 2022-07-10         4.0
    4     Amy 2022-07-11         5.0
    ..    ...        ...         ...
    47     Gu 2022-07-07         2.0
    48     Gu 2022-07-11         5.0
    49     Gu 2022-07-12         2.0
    50     Gu 2022-07-13         4.0
    51     Gu 2022-07-14         3.0
    <BLANKLINE>
    [52 rows x 3 columns]

The solution can be transformed into alternative output formats using standard
pandas operations. For example, the shift assignments could be pivoted to
produce a wide-format table displaying a readable roster. Alternatively, one
could use pandas I/O functions to push the solution to another system or service
for further processing.

.. doctest:: workforce
    :options: +NORMALIZE_WHITESPACE

    >>> shifts_table = pd.pivot_table(
    ...     assigned_shifts.assign(value=1),
    ...     values="value",
    ...     index="Shift",
    ...     columns="Worker",
    ...     fill_value="-",
    ... ).replace({1.0: "Y"})
    >>> with pd.option_context('display.max_rows', 15):
    ...     print(shifts_table)
    Worker     Amy Bob Cathy Dan Ed Fred Gu
    Shift
    2022-07-01   -   Y     -   -  -    Y  Y
    2022-07-02   -   -     -   Y  Y    -  -
    2022-07-03   Y   -     Y   Y  Y    -  -
    2022-07-04   -   -     Y   -  Y    -  -
    2022-07-05   Y   Y     Y   Y  -    -  Y
    2022-07-06   -   Y     -   Y  -    Y  Y
    2022-07-07   Y   -     Y   -  Y    -  Y
    2022-07-08   -   Y     -   -  Y    -  -
    2022-07-09   -   -     -   -  Y    Y  -
    2022-07-10   Y   -     Y   Y  -    -  -
    2022-07-11   Y   Y     -   Y  -    -  Y
    2022-07-12   Y   -     Y   Y  -    Y  Y
    2022-07-13   Y   Y     Y   Y  Y    Y  Y
    2022-07-14   Y   -     Y   -  Y    Y  Y

Enforcing Breaks
----------------

The approach above is likely too simple for longer rosters, since the number of
shifts assigned to each worker is only constrained over the entire time period
of the roster. Realistically, this requirement may need to be enforced on a
rolling basis, for example a worker may only be allowed to be assigned four
shifts in any given five day period (i.e. one rostered-off day). This is
enforced using the ``limit_window`` keyword argument. If this optional
argument is provided, the ``worker_limits`` constraint will be enforced over
every rolling window of the given time period, instead of over the entire roster.

.. doctest:: workforce
    :options: +NORMALIZE_WHITESPACE +ELLIPSIS

    >>> worker_limits = pd.DataFrame(dict(
    ...     Worker=data.worker_limits["Worker"],
    ...     MinShifts=0,
    ...     MaxShifts=4,
    ... ))
    >>> worker_limits
      Worker  MinShifts  MaxShifts
    0    Amy          0          4
    1    Bob          0          4
    2  Cathy          0          4
    3    Dan          0          4
    4     Ed          0          4
    5   Fred          0          4
    6     Gu          0          4
    >>> assigned_shifts = solve_workforce_scheduling(
    ...     preferences=data.preferences,
    ...     shift_requirements=data.shift_requirements,
    ...     worker_limits=worker_limits,
    ...     limit_window=pd.Timedelta("5D"),
    ...     silent=True,
    ... )
    >>> shifts_table = pd.pivot_table(
    ...     assigned_shifts.assign(value=1),
    ...     values="value",
    ...     index="Shift",
    ...     columns="Worker",
    ...     fill_value="-",
    ... ).replace({1.0: "Y"})
    >>> with pd.option_context('display.max_rows', 15):
    ...     print(shifts_table)
    Worker     Amy Bob Cathy Dan Ed Fred Gu
    Shift
    2022-07-01   -   Y     -   -  -    Y  Y
    2022-07-02   -   -     -   Y  Y    -  -
    2022-07-03   -   -     Y   Y  Y    -  Y
    2022-07-04   -   -     Y   -  Y    -  -
    2022-07-05   Y   -     Y   Y  Y    -  Y
    2022-07-06   -   Y     -   Y  -    Y  Y
    2022-07-07   Y   -     Y   -  Y    -  Y
    2022-07-08   -   -     -   -  Y    -  Y
    2022-07-09   -   -     -   -  Y    Y  -
    2022-07-10   Y   -     Y   Y  -    -  -
    2022-07-11   -   Y     -   Y  Y    -  Y
    2022-07-12   Y   -     Y   Y  -    Y  Y
    2022-07-13   Y   Y     Y   Y  Y    Y  Y
    2022-07-14   Y   -     Y   -  Y    Y  Y

Observe that the schedule has been adjusted to avoid any worker working
more than 5 consecutive days.
