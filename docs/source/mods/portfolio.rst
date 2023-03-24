Portfolio Optimization
======================

Portfolio optimization is concerned with allocation of wealth into assets.


Problem Specification
---------------------

Our methods use risk and return estimators.

.. tabs::

    .. tab:: Domain-Specific Description

        We consider a multi-periods portfolio optimization problem where ...

    .. tab:: Optimization Model

        Give the mathematical programming formulation of the problem here.
        The risk is expressed as

        .. math::

            x^\top \Sigma x

            x^\mathrm{T} \Sigma x


The risk estimator is given through a time series.

.. testsetup:: mod

    # Set pandas options for displaying dataframes, if needed
    import pandas as pd
    pd.options.display.max_rows = 10
    pd.options.display.max_columns = 7

.. tabs::

    .. tab:: ``stocks``

        For a number of stocks we have historic values.

        .. doctest:: mod
            :options: +NORMALIZE_WHITESPACE

            >>> from gurobi_optimods import datasets
            >>> data = datasets.load_portfolio()
            >>> data
                       AA        BB        CC  ...        HH        II        JJ
            0   -0.000601  0.002353 -0.075234  ...  0.060737 -0.012869 -0.022137
            1    0.019177  0.050008  0.041345  ... -0.026674  0.009876  0.012809
            2   -0.020333  0.026638 -0.038999  ... -0.023153 -0.007007 -0.038034
            3    0.001421 -0.053813 -0.013347  ... -0.012348  0.018736 -0.058373
            4    0.008648  0.114836  0.003617  ...  0.064090 -0.011153  0.024333
            ..        ...       ...       ...  ...       ...       ...       ...
            240  0.007382 -0.000724 -0.002444  ... -0.007654  0.018015 -0.017135
            241 -0.003362 -0.106913 -0.082365  ...  0.040787 -0.023020 -0.067792
            242  0.004359 -0.029763 -0.041986  ...  0.014522 -0.017109 -0.060247
            243 -0.018402 -0.054211 -0.075788  ... -0.013557  0.022576 -0.036793
            244 -0.016237  0.015580 -0.026970  ... -0.005893 -0.013456 -0.032203
            <BLANKLINE>
            [245 rows x 10 columns]

        In the model, this corresponds to ...


Code
----

Self contained code example to run the mod from an example dataset. Example
datasets should bd included in the ``gurobi_optimods.datasets`` module for
easy access by users.

.. testcode:: mod

    import pandas as pd

    from gurobi_optimods.datasets import load_portfolio
    from gurobi_optimods.portfolio import solve_mod


    data = load_portfolio()
    solution = solve_mod(data)

..  A snippet of the Gurobi log output here won't show in the rendered page,
    but serves as a doctest to make sure the code example runs. The ... lines
    are meaningful here, they will match anything in the output test.

.. testoutput:: mod
    :hide:

    stocks


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
