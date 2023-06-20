Workforce Scheduling
====================

Workforce scheduling is an extremely widely-used application of optimization in
practice. It involves balancing many competing concerns such as worker
availability, cost, and preferences, shift coverage requirements, conditions on
consecutive shifts or rest breaks, etc. Implementation can become quite
involved as worker requirements and entitlements become more complex.

This Mod implements several basic variants of workforce scheduling. The initial
example is deliberately simple, while later sections progressively add more
complexity and enforce additional requirements on the generated schedule. If the
initial example appears too simple for your use case, please, read on!

Problem Specification
---------------------

This first example covers a simple case for a business developing a two-week
roster. Each shift requires a given number of workers who have identical skills
and productivity. In other words, any work can cover any shift, and all workers
are considered equivalent. Workers provide their availability for shifts.
Rest requirements and minimum work entitlements are handled through upper and
lower limits, respectively, on the number of shifts a worker is rostered on
in the schedule. Optionally, preferences can be provided, in which case the
scheduler aims to maximize satisfaction by finding a feasible roster that
maximizes the sum of preference scores of the assigned shifts.

The workforce scheduling Mod takes the following three pandas dataframes as input:

* The ``availability`` dataframe has two columns: ``Worker`` and
  ``Shift``. A row in this dataframe indicates that the given worker is
  available to work the given shift. If the optional ``preferences`` argument
  is provided, it refers to an additional column in the ``availability``
  dataframe that provides numerical preference data.
* The ``shift_requirements`` dataframe has two columns: ``Shift`` and
  ``Required``. Each row specifies the number of workers required for a given
  shift. There must be one row for every unique shift mentioned in
  ``availability["Shift"]``.
* The ``worker_limits`` dataframe has three columns: ``Worker``,
  ``MinShifts``, and ``MaxShifts``. Each row specifies the minimum and maximum
  number of shifts the given worker may be assigned in the schedule. There
  must be one row for every unique worker mentioned in ``availability["Worker"]``.

When ``solve_workforce_scheduling`` is called, a model is formulated and solved
immediately using Gurobi. Workers will be assigned only to shifts they are
available for, in such a way that all requirements are covered, minimum and
maximum shift numbers are respected, and the total sum of worker preference
scores is maximised. If ``preferences=None``, preferences are not considered and any
feasible schedule will be returned.

The returned assignment dataframe is a subset of the availability dataframe,
with the same columns. A row in the returned dataframe specifies that the given
worker has been assigned to the given shift.

.. dropdown:: Background: Mathematical Model

    A set of shifts :math:`S` is to be covered using a set of workers
    :math:`W`. Workers :math:`w \in W_{s} \subseteq W` are available to work
    a given shift `s`, and have a preference :math:`p_{ws}` for each
    assigned shift. Shift :math:`s` requires :math:`r_{s}` workers assigned.
    Each worker must be assigned between :math:`l_{w}` and :math:`u_{w}`
    shifts in total. The model is defined on binary variables :math:`x_{ws}`
    that satisfy the condition

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

A Simple Example
----------------

This section shows the simplest possible input dataset for the Mod, comprising
seven workers covering daily shifts over a two week period, each with defined
availability. All input data is given as pandas dataframes, following the layout
described in the :ref:`mods/workforce:problem specification` above. The tabs
below show example data for each frame.

.. testsetup:: workforce

    # Set pandas options
    import pandas as pd
    pd.options.display.max_rows = 10

.. tabs::

    .. tab:: ``availability``

        The following example table lists worker availability and preferences.
        For example, Siva is available on May 2nd, 3rd, 5th, and so on, with a
        stronger preference to be assigned the shift on the 5th. To use the
        preference data, the optional argument ``preferences="Preference"`` must
        be supplied.

        .. doctest:: workforce
            :options: +NORMALIZE_WHITESPACE

            >>> from gurobi_optimods import datasets
            >>> data = datasets.load_workforce()
            >>> data.availability
                 Worker      Shift  Preference
            0      Siva 2023-05-02         2.0
            1      Siva 2023-05-03         2.0
            2      Siva 2023-05-05         5.0
            3      Siva 2023-05-07         3.0
            4      Siva 2023-05-09         2.0
            ..      ...        ...         ...
            67  Pauline 2023-05-10         4.0
            68  Pauline 2023-05-11         5.0
            69  Pauline 2023-05-12         2.0
            70  Pauline 2023-05-13         4.0
            71  Pauline 2023-05-14         3.0
            <BLANKLINE>
            [72 rows x 3 columns]

        In the mathematical model, the worker-shift pairings enumerate all
        possible members of the set
        :math:`\lbrace (w, s) \mid s \in S, w \in W_s \rbrace`, and the
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
            0  2023-05-01         3
            1  2023-05-02         2
            2  2023-05-03         4
            3  2023-05-04         2
            4  2023-05-05         5
            ..        ...       ...
            9  2023-05-10         3
            10 2023-05-11         4
            11 2023-05-12         5
            12 2023-05-13         7
            13 2023-05-14         5
            <BLANKLINE>
            [14 rows x 2 columns]

        In the mathematical model, this table provides the values :math:`r_s`.

    .. tab:: ``worker_limits``

        The following example table lists the minimum and maximum number of
        shifts in the planning period that each worker is entitled to.

        .. doctest:: workforce
            :options: +NORMALIZE_WHITESPACE

            >>> from gurobi_optimods import datasets
            >>> data = datasets.load_workforce()
            >>> data.worker_limits
                Worker  MinShifts  MaxShifts
            0     Siva          6          8
            1  Ziqiang          6          7
            2  Matsumi          6          8
            3    Femke          5          8
            4  Vincent          6          8
            5   Marisa          5          8
            6  Pauline          6          8

        In the mathematical model, this table provides the values :math:`l_w`
        and :math:`u_w`.

The example code below solves the workforce scheduling problem for the above
dataset. The dataset can be imported from the ``gurobi-optimods.datasets`` module.

.. testcode:: workforce

    from gurobi_optimods.datasets import load_workforce
    from gurobi_optimods.workforce import solve_workforce_scheduling

    # Load example data
    data = load_workforce()

    # Solve the Mod, get back a schedule
    assigned_shifts = solve_workforce_scheduling(
        availability=data.availability,
        shift_requirements=data.shift_requirements,
        worker_limits=data.worker_limits,
        preferences="Preference",
    )

.. testoutput:: workforce
    :hide:

    ...
    Optimize a model with 28 rows, 72 columns and 216 nonzeros
    ...
    Best objective 1.850000000000e+02, best bound 1.850000000000e+02, gap 0.0000%

Inspecting the Solution
-----------------------

The solution to this workforce scheduling problem is a selection of shift
assignments. The returned dataframe is a subset of the original ``availability``
dataframe.

.. doctest:: workforce
    :options: +NORMALIZE_WHITESPACE

    >>> assigned_shifts
          Worker      Shift  Preference
    0       Siva 2023-05-03         2.0
    1       Siva 2023-05-05         5.0
    2       Siva 2023-05-07         3.0
    3       Siva 2023-05-10         4.0
    4       Siva 2023-05-11         5.0
    ..       ...        ...         ...
    47   Pauline 2023-05-07         2.0
    48   Pauline 2023-05-11         5.0
    49   Pauline 2023-05-12         2.0
    50   Pauline 2023-05-13         4.0
    51   Pauline 2023-05-14         3.0
    <BLANKLINE>
    [52 rows x 3 columns]

The solution can be transformed into alternative output formats using standard
pandas operations. For example, the shift assignments could be pivoted to
produce a wide-format table displaying a readable roster. Alternatively, one
could use pandas I/O functions to push the solution to another system or service
for further processing.

.. doctest:: workforce
    :options: +NORMALIZE_WHITESPACE

    >>> import pandas as pd
    >>> shifts_table = pd.pivot_table(
    ...     assigned_shifts.assign(value=1),
    ...     values="value",
    ...     index="Shift",
    ...     columns="Worker",
    ...     fill_value="-",
    ... ).replace({1.0: "Y"})
    >>> with pd.option_context('display.max_rows', 15):
    ...     print(shifts_table)
    Worker     Femke Marisa Matsumi Pauline Siva Vincent Ziqiang
    Shift
    2023-05-01     -      Y       -       Y    -       -       Y
    2023-05-02     Y      -       -       -    -       Y       -
    2023-05-03     Y      -       Y       -    Y       Y       -
    2023-05-04     -      -       Y       -    -       Y       -
    2023-05-05     Y      -       Y       Y    Y       -       Y
    2023-05-06     Y      Y       -       Y    -       -       Y
    2023-05-07     -      -       Y       Y    Y       Y       -
    2023-05-08     -      -       -       -    -       Y       Y
    2023-05-09     -      Y       -       -    -       Y       -
    2023-05-10     Y      -       Y       -    Y       -       -
    2023-05-11     Y      -       -       Y    Y       -       Y
    2023-05-12     Y      Y       Y       Y    Y       -       -
    2023-05-13     Y      Y       Y       Y    Y       Y       Y
    2023-05-14     -      Y       Y       Y    Y       Y       -

Enforcing Breaks
----------------

The approach above is likely too simple for longer rosters, since the number of
shifts assigned to each worker is only constrained over the entire time period
of the roster. Realistically, this requirement may need to be enforced on a
rolling basis.  For example, a worker may only be allowed to be assigned four
shifts in any given five day period (i.e. one rostered-off day). This is
enforced using the ``limit_window`` keyword argument. If this optional
argument is provided, the ``worker_limits`` constraint will be enforced over
a rolling window of the given duration, rather than over the entire roster
duration.

.. doctest:: workforce
    :options: +NORMALIZE_WHITESPACE +ELLIPSIS

    >>> worker_limits = pd.DataFrame(dict(
    ...     Worker=data.worker_limits["Worker"],
    ...     Window=pd.Timedelta("5D"),
    ...     MinShifts=0,
    ...     MaxShifts=4,
    ... ))
    >>> worker_limits
        Worker Window  MinShifts  MaxShifts
    0     Siva 5 days          0          4
    1  Ziqiang 5 days          0          4
    2  Matsumi 5 days          0          4
    3    Femke 5 days          0          4
    4  Vincent 5 days          0          4
    5   Marisa 5 days          0          4
    6  Pauline 5 days          0          4

The above data specifies that all workers have identical requirements to work at
most four shifts in any given five day period, with no minimum number of shifts
required. When solving this variant of the problem, ``rolling_limits`` must be
set to ``True`` to enforce the new requirement.

.. doctest:: workforce
    :options: +NORMALIZE_WHITESPACE +ELLIPSIS

    >>> assigned_shifts = solve_workforce_scheduling(
    ...     availability=data.availability,
    ...     shift_requirements=data.shift_requirements,
    ...     worker_limits=worker_limits,
    ...     preferences="Preference",
    ...     rolling_limits=True,
    ...     verbose=False,
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
    Worker     Femke Marisa Matsumi Pauline Siva Vincent Ziqiang
    Shift
    2023-05-01     -      Y       -       Y    -       -       Y
    2023-05-02     Y      -       -       -    -       Y       -
    2023-05-03     Y      -       Y       Y    -       Y       -
    2023-05-04     -      -       Y       -    -       Y       -
    2023-05-05     Y      -       Y       Y    Y       Y       -
    2023-05-06     Y      Y       -       Y    -       -       Y
    2023-05-07     -      -       Y       Y    Y       Y       -
    2023-05-08     -      -       -       Y    -       Y       -
    2023-05-09     -      Y       -       -    -       Y       -
    2023-05-10     Y      -       Y       -    Y       -       -
    2023-05-11     Y      -       -       Y    -       Y       Y
    2023-05-12     Y      Y       Y       Y    Y       -       -
    2023-05-13     Y      Y       Y       Y    Y       Y       Y
    2023-05-14     -      Y       Y       Y    Y       Y       -

Notice that Siva's shifts have been adjusted so as to avoid any worker working
more than five consecutive days.

Further Requirements
--------------------

As mentioned in the introduction, this Mod implements some basic cases of
workforce scheduling, and is limited in scope. However, similar modelling
approaches to those described here can be applied to handle more complex
requirements. For further information, see :footcite:t:`ERNST20043` (among many,
many other references on the topic).

.. footbibliography::
