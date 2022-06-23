Workforce Scheduling
====================

A little background on the proud history of mathprog in workforce scheduling.

Also data science.


Problem Specification
---------------------

General problem description. Tabs alternate between describing the input dataframe structures and defining the math model.

.. tabs::

    .. tab:: Description

        * Each row in availability dataframe specifies that an employee is available to work a given shift
        * Pay dataframe specifies the cost of having employee x work a shift
        * Requirements dataframe specifies the number of workers required on each shift
        * Employees will be assigned only to shifts they are available for, in such a way that all requirements are covered while total cost of covering all shifts is minimised
        * Each row in the returned assignment dataframe specifies that an employee has been assigned a given shift

    .. tab:: Maths

        I should define some terms here...

        .. math::

            \min \sum_i pay_{shift} x_{shift}

            \sum_i x_{shift} = req_{shift}

|

Code
----

Alternate between the code required to run the model from the store vs how to implement directly in gurobipy. If you use nupstup, all the gurobi internals are handled for you; users interact with the 'solver' by passing dataframes to a given spec and receiving a dataframe as output. If you instead peek under the hood and use gurobipy, you have more options to extend the model with additional constraints and data.

.. tabs::
    .. tab:: Nupstup

        .. literalinclude:: ../../examples/workforce_nupstup.py
            :linenos:
            :caption: workforce_nupstup.py

    .. tab:: Gurobipy

        .. literalinclude:: ../../examples/workforce_gurobipy.py
            :linenos:
            :caption: workforce_gurobipy.py


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
