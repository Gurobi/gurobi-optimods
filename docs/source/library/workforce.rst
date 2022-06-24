Workforce Scheduling
====================

A little background on the proud history of mathprog in workforce scheduling.

Also data science.

Problem Specification
---------------------

Consider a service business, like a restaurant, that develops its workforce plans for the next two weeks (considering a 7-day week). The service requires only one set of skills. There are a number of employed workers with the same set of skills and with identical productivity that are available to work on some of the days during the two-week planning horizon. There is only one shift per workday. Each shift may have different resource (worker) requirements on each workday. The service business wants to minimize the cost of covering all shifts with available workers, based on their rates of pay and availability.

.. tabs::

    .. tab:: Data Specification

        The workforce scheduling model takes the following inputs:

        * The ``availability`` dataframe has two columns: ``Worker`` and ``Shift``. Each row in dataframe specifies that the given worker is available to work the given shift.
        * The ``shift_requirements`` series is indexed by shift and specifies the number of workers required for that shift. There should be one entry for every unique worker in ``availability["Workers"]``.
        * The ``pay_rates`` series is indexed by worker and specifies the pay per shift worked for that worker. There should be one entry for every unique shift in ``availability["Shift"]``.

        When ``solve_workforce_scheduling`` is called, a model is formulated and solved immediately using Gurobi. Workers will be assigned only to shifts they are available for, in such a way that all requirements are covered while total cost of covering all shifts is minimised.

        The returned assignment dataframe is a subset of the availability dataframe, with the same columns. Each row specifies that the given worker has been assigned the given shift.

    .. tab:: Mathematical Model

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

.. tabs::

    .. tab:: ``availability``

        Amy is available for a shift on Tuesday 2nd, etc, etc

        .. code::

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

            [72 rows x 2 columns]

        In the mathematical model, this models the set :math:`\lbrace (w, s) \mid s \in S, w \in W_s \rbrace`.

    .. tab:: ``shift_requirements``

        Shift on Monday 1st requires 3 workers, etc, etc

        .. code::

                Required      Shift
            0          3 2022-07-01
            1          2 2022-07-02
            2          4 2022-07-03
            3          2 2022-07-04
            4          5 2022-07-05
            5          4 2022-07-06
            6          4 2022-07-07
            7          2 2022-07-08
            8          2 2022-07-09
            9          3 2022-07-10
            10         4 2022-07-11
            11         5 2022-07-12
            12         7 2022-07-13
            13         5 2022-07-14

        In the mathematical model, this models the values :math:`r_s`.

    .. tab:: ``pay_rates``

        Bob is the most expensive worker ...

        .. code::

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

Alternate between the code required to run the model from the store vs how to implement directly in gurobipy. If you use nupstup, all the gurobi internals are handled for you; users interact with the 'solver' by passing dataframes to a given spec and receiving a dataframe as output. If you instead peek under the hood and use gurobipy, you have more options to extend the model with additional constraints and data.

.. tabs::
    .. tab:: Nupstup

        .. literalinclude:: ../../examples/workforce_nupstup.py
            :linenos:

    .. tab:: Gurobipy

        .. literalinclude:: ../../examples/workforce_gurobipy.py
            :linenos:


Both codes construct the same model and give the same result. The model is solved as a linear program by Gurobi.

.. collapse:: View Gurobi logs

    .. literalinclude:: ../../examples/workforce.log
        :language: none

|

Solution
--------

Solution is a selection of shift assignments.

.. code-cell:: python
    :execution-count: 1

    assigned_shifts

.. output-cell::
    :execution-count: 1

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

Use pandas functions to create a shift allocation table for added prettiness.

.. code-cell:: python
    :execution-count: 2

    shifts_table = pd.pivot_table(
        assigned_shifts.assign(value=1),
        values="value",
        index="Shift",
        columns="Worker",
        fill_value="-",
    ).replace({1.0: "Y"})

    shifts_table

.. output-cell::
    :execution-count: 2

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
