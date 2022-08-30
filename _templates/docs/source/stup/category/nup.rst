My New Nup
==========

A little background on the proud history of mathprog in this field.

Also data science.

Problem Specification
---------------------

Give a brief overview of the problem being solved.

.. tabs::

    .. tab:: Domain-Specific Description

        Give a definition of the problem in the language of the domain expert.

    .. tab:: Mathematical Model

        Give the mathematical programming formulation of the problem here.

Give examples of the various input data structures, if appropriate. These outputs should be fixed, so use doctests where possible.

.. testsetup:: nup

    # Set pandas options for displaying dataframes, if needed
    import pandas as pd
    pd.options.display.max_rows = 10

.. tabs::

    .. tab:: ``availability``

        Give interpretation of input data.

        .. doctest:: nup
            :options: +NORMALIZE_WHITESPACE

            >>> pd.read_feather("examples/data/availability.feather")
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

        In the mathematical model, this models ...

    .. tab:: ``shift_requirements``

|

Code
----

Tab between the code required to run the model from the store vs how to implement directly in gurobipy. If you use nupstup, all the gurobi internals are handled for you; users interact with the 'solver' by passing dataframes to a given spec and receiving a dataframe as output. If you instead peek under the hood and use gurobipy, you have more options to extend the model with additional constraints and data.

.. These paths need to be changed to point to your example scripts

.. tabs::
    .. tab:: nupstup function

        .. literalinclude:: ../../../examples/nup/nupstup.py
            :linenos:

    .. tab:: gurobipy model

        .. literalinclude:: ../../../examples/nup/gurobipy.py
            :linenos:


Both codes construct the same model and give the same result. The model is solved as a LP/MIP/QP/etc by Gurobi.

.. collapse:: View Gurobi logs

    .. code-block:: text

        Gurobi Optimizer version 9.5.1 build v9.5.1rc2
        ...

|

Solution
--------

Show the solution. Use doctests if possible (i.e. the solution must be stable enough). Otherwise, just display it somehow.

.. This import line needs to be changed to import any results you need

.. testcode:: nup
    :hide:

    from examples.<nup>.nupstup import some_result

.. testoutput:: nup
    :hide:

    Gurobi Optimizer version 9.5.1 build v9.5.1rc2
    ...

.. doctest:: nup
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
