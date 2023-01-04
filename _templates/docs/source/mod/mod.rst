My New Mod
==========

A little background on the proud history of mathprog in this field.

Also data science.

Problem Specification
---------------------

Give a brief overview of the problem being solved.

.. tabs::

    .. tab:: Domain-Specific Description

        Give a definition of the problem in the language of the domain expert.

    .. tab:: Optimization Model

        Give the mathematical programming formulation of the problem here.

Give examples of the various input data structures. These inputs should be fixed,
so use doctests where possible.

.. testsetup:: mod

    # Set pandas options for displaying dataframes, if needed
    import pandas as pd
    pd.options.display.max_rows = 10

.. tabs::

    .. tab:: ``availability``

        Give interpretation of input data.

        .. doctest:: mod
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

        In the model, this corresponds to ...

    .. tab:: ``shift_requirements``

        Another bit of input data (perhaps a secondary table)

|

Code
----

Self contained code example to run the mod from an example dataset. Example
datasets should bd included in the ``gurobi_optimods.datasets`` module for
easy access by users.

.. testcode:: mod

    import pandas as pd

    from gurobi_optimods.datasets import load_mod_data
    from gurobi_optimods.mod import solve_mod


    data = load_mod_data()
    solution = solve_mod(data.table1, data.table2)

..  A snippet of the Gurobi log output here won't show in the rendered page,
    but serves as a doctest to make sure the code example runs. The ... lines
    are meaningful here, they will match anything in the output test.

.. testoutput:: mod
    :hide:

    ...
    Optimize a model with 14 rows, 72 columns, and 72 nonzeros
    ...
    Optimal objective
    ...

The model is solved as an LP/MIP/QP by Gurobi.

..  You can include the full Gurobi log output here for the curious reader.
    It will be visible as a collapsible section.

.. collapse:: View Gurobi Logs

    .. code-block:: text

        Gurobi Optimizer version 9.5.1 build v9.5.1rc2 (mac64[x86])
        Optimize a model with ...
        Best obj ... Best bound ...

|

Solution
--------

Show the solution. One way is to use doctests to display simple shell outputs
(see the workforce example). Another option is to include and display figures
(see the graph matching examples).
