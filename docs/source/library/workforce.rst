Workforce Scheduling
====================

Problem Specification
---------------------

General problem description. Tabs alternate between describing the input
dataframe structures and defining the math model.

.. tabs::

    .. group-tab:: Description

        each row in availability dataframe specifies that an employee is
        available to work a given shift ...

    .. group-tab:: Maths

        .. math::

            \sum_i x_{shift} = req_{shift}


Code
----

Alternate between the code required to run the model from the store
vs how to implement directly in gurobipy (needed when a user wants
to extend).

.. tabs::
    .. group-tab:: Nupstup

        .. literalinclude:: ../../examples/workforce_nupstup.py
            :linenos:
            :caption: workforce_nupstup.py

    .. group-tab:: Gurobipy

        .. literalinclude:: ../../examples/workforce_gurobipy.py
            :linenos:
            :caption: workforce_gurobipy.py


Logging output, model handled as a pure LP.
Maybe make this an expandable thing.

.. literalinclude:: ../../examples/workforce.log
    :caption: Gurobi logging output


Solution is a selection of shift assignments.

.. literalinclude:: ../../examples/workforce.out
    :caption: Solution
