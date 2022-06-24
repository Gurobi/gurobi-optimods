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

        * The ``availability`` dataframe has two columns: ``Workers`` and ``Shift``. Each row in dataframe specifies that the given worker is available to work the given shift.
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

               Workers  Shift
            0      Amy   Tue2
            1      Amy   Wed3
            2      Amy   Fri5
            3      Amy   Sun7
            4      Amy   Tue9
            ..     ...    ...
            67      Gu  Wed10
            68      Gu  Thu11
            69      Gu  Fri12
            70      Gu  Sat13
            71      Gu  Sun14

            [72 rows x 2 columns]

    .. tab:: ``shift_requirements``

        Shift on Monday 1st requires 3 workers, etc, etc

        .. code::

            Shift
            Mon1     3
            Tue2     2
            Wed3     4
            Thu4     2
            Fri5     5
            ..      ...
            Wed10    3
            Thu11    4
            Fri12    5
            Sat13    7
            Sun14    5
            Name: Req, dtype: int64

    .. tab:: ``pay_rates``

        Bob is the most expensive worker ...

        .. code::

            Workers
            Amy      10
            Bob      12
            Cathy    10
            Dan       8
            Ed        8
            Fred      9
            Gu       11
            Name: Pay, dtype: int64

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

       Workers  Shift
    1      Amy   Wed3
    2      Amy   Fri5
    3      Amy   Sun7
    ...
    64      Gu   Sun7
    69      Gu  Fri12
    70      Gu  Sat13

Use pandas functions to create a shift allocation table for added prettiness.

.. code-cell:: python
    :execution-count: 2

    shifts_table = pd.pivot_table(
        assigned_shifts.assign(value=1),
        values="value",
        index="Shift",
        columns="Workers",
        fill_value="-",
    ).replace({1.0: "Y"})

    shifts_table

.. output-cell::
    :execution-count: 2

    Workers Amy Bob Cathy Dan Ed Fred Gu
    Shift
    Fri12     Y   -     Y   Y  -    Y  Y
    Fri5      Y   -     Y   Y  Y    -  Y
    Mon1      -   -     -   -  Y    Y  Y
    Mon8      -   -     -   Y  Y    -  -
    Sat13     Y   Y     Y   Y  Y    Y  Y
    Sat6      -   Y     -   Y  -    Y  Y
    Sun14     Y   -     Y   Y  Y    Y  -
    Sun7      Y   -     Y   -  Y    -  Y
    Thu11     Y   -     Y   Y  Y    -  -
    Thu4      -   -     Y   -  Y    -  -
    Tue2      -   -     -   Y  Y    -  -
    Tue9      -   -     -   Y  Y    -  -
    Wed10     Y   -     Y   Y  -    -  -
    Wed3      Y   -     -   Y  Y    Y  -
