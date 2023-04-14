ACOPF
=====

Prior to specifying the mathematical structure of the problem it is important to present basic concepts and terminology used by the power systems community.

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

    from gurobi_optimods.opf import solve_opf_model, read_settings_from_file, read_case_from_file
    from gurobi_optimods.datasets import load_caseopf

    settings = {"doac": True, "use_ef": True}
    # load path to case file
    casefile = load_caseopf("9")
    # read case file and return a case dictionary
    case = read_case_from_file(casefile)
    # solve opf model and return a solution and the final objective value
    solution, objval = solve_opf_model(settings, case)

..  A snippet of the Gurobi log output here won't show in the rendered page,
    but serves as a doctest to make sure the code example runs. The ... lines
    are meaningful here, they will match anything in the output test.

.. testoutput:: mod
    :hide:

    ...
    Optimize a model with 73 rows, 107 columns and 208 nonzeros
    ...
    Optimal solution found (tolerance 1.00e-03)
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
(see the workforce example). This can be done simply by pasting outputs
directly from a python shell. Another option is to include and display figures
(see the graph matching examples).

.. doctest:: mod
    :options: +NORMALIZE_WHITESPACE

    >>>
