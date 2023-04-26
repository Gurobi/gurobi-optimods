Optimal mean-variance portfolios
================================

Portfolio optimization is concerned with the allocation of wealth into assets (such as stocks, bonds, commodities, etc.). This mod returns portfolios on the efficient frontier given expected returns and variances.


Problem Specification
---------------------

Our methods use risk and return estimators.

.. tabs::

    .. tab:: Domain-Specific Description

        We consider a single-period portfolio optimization problem where want to allocate wealth into :math:`n` risky assets. The returned portfolio :math:`x` is an efficient mean-variance portfolio given returns :math:`\mu`, covariance :math:`\Sigma` and risk aversion :math:`\gamma`.


    .. tab:: Mathematical Description

        .. math::

            \begin{alignat}{2}
            \max \quad        & \mu^\top x - \tfrac12 \gamma\ x^\top\Sigma x \\
            \mbox{s.t.} \quad & 1^\top x = 1 \\
            \end{alignat}

        * :math:`\mu` is the vector of expected returns.
        * :math:`\Sigma` is the return covariance matrix.
        * :math:`x` is the portfolio where :math:`x_i` denotes the fraction of wealth invested in the risky asset :math:`i`.
        * :math:`\gamma\geq0` is the risk aversion coefficient.


This description refers only to the simple base model.  Further down in Section
`Enforcing more portfolio features`_ we explain how to enforce additional
fatures, such as leverage or transaction fees.


Data Specification
------------------

The mean-variance portfolio optimization model takes the following inputs:

* The covariance matrix :math:`\Sigma` can be given as a pandas Dataframe or a numpy array.
* The return estimator :math:`\mu` can be given as a pandas Series or a numpy array.


The returned allocation vector :math:`x` is a pandas Series or a numpy ndarray (depending on the input types).


Example
-------

Here we derive the matrix :math:`\Sigma` and the vector :math:`\mu` from a time series of logarithmic historic returns:

.. testsetup:: mod

    # Set pandas options for displaying dataframes, if needed
    import pandas as pd
    pd.options.display.max_rows = 10
    pd.options.display.max_columns = 7

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

The columns of this dataframe represent the individual assets ("AA", "BB", ..., "JJ") while the rows represent the historic time steps. We use pandas functionality to compute a simple mean estimator and corresponding covariance from this dataframe:


.. testcode:: mod

    import pandas as pd

    from gurobi_optimods.datasets import load_portfolio
    from gurobi_optimods.portfolio import MeanVariancePortfolio

    data = load_portfolio()
    Sigma = data.cov()
    mu = data.mean()
    gamma = 100.0

    mvp = MeanVariancePortfolio(Sigma, mu)
    x = mvp.efficient_portfolio(gamma)

.. testoutput:: mod
    :hide:

    ...
    Optimize a model with 42 rows, 50 columns and 100 nonzeros
    ...
    Model has 55 quadratic objective terms
    ...



..  You can include the full Gurobi log output here for the curious reader.
    It will be visible as a collapsible section.

.. collapse:: View Gurobi Logs

    .. code-block:: text

        Gurobi Optimizer version 10.0.1 build v10.0.1rc0 (mac64[rosetta2])

        CPU model: Apple M1
        Thread count: 8 physical cores, 8 logical processors, using up to 8 threads

        Optimize a model with 1 rows, 10 columns and 10 nonzeros
        Model fingerprint: 0x7edd9de0
        Model has 55 quadratic objective terms
        Coefficient statistics:
        Matrix range     [1e+00, 1e+00]
        Objective range  [7e-04, 1e-02]
        QObjective range [7e-06, 2e-03]
        Bounds range     [0e+00, 0e+00]
        RHS range        [1e+00, 1e+00]
        Presolve time: 0.01s
        Presolved: 1 rows, 10 columns, 10 nonzeros
        Presolved model has 55 quadratic objective terms
        Ordering time: 0.00s

        Barrier statistics:
        Free vars  : 9
        AA' NZ     : 4.500e+01
        Factor NZ  : 5.500e+01
        Factor Ops : 3.850e+02 (less than 1 second per iteration)
        Threads    : 1

                          Objective                Residual
        Iter       Primal          Dual         Primal    Dual     Compl     Time
           0  -2.08348238e+05  2.08383773e+05  1.00e+04 1.43e-02  1.00e+06     0s
           1  -1.91482256e-01  4.99463850e+02  1.08e+01 9.88e-09  1.12e+03     0s
           2  -1.94725618e-02  4.56374984e+02  1.08e-05 9.88e-15  4.56e+01     0s
           3  -1.94685319e-02  4.71448851e-01  8.14e-10 1.39e-17  4.91e-02     0s
           4  -1.63767350e-02  1.14105476e-02  2.04e-11 6.94e-18  2.78e-03     0s
           5  -7.58352892e-03  1.59186002e-04  3.89e-16 2.08e-17  7.74e-04     0s
           6  -5.59221914e-03 -4.72740622e-03  1.67e-16 6.94e-18  8.65e-05     0s
           7  -5.18009820e-03 -5.10195350e-03  9.30e-16 1.04e-17  7.81e-06     0s
           8  -5.12692872e-03 -5.12414839e-03  6.11e-16 3.47e-18  2.78e-07     0s
           9  -5.12425311e-03 -5.12424841e-03  4.84e-15 6.94e-18  4.70e-10     0s

        Barrier solved model in 9 iterations and 0.00 seconds (0.00 work units)
        Optimal objective -5.12425311e-03


Solution
--------

The returned Series contains the relative investment for each asset;
here the solution suggests to spread the investments over five positions
(AA, DD, GG, HH, II).  The other allocations are negligible.

.. doctest:: mod
    :options: +NORMALIZE_WHITESPACE

    >>> x
    AA    4.236507e-01
    BB    1.743570e-07
    CC    7.573610e-10
    DD    2.430104e-01
    EE    1.017732e-07
    FF    2.760531e-09
    GG    2.937307e-02
    HH    2.350833e-01
    II    6.888222e-02
    JJ    1.248442e-08
    dtype: float64


Enforcing more portfolio features
---------------------------------

A number of additional restrictions can be placed on the returned optimal
portfolio, such as transaction fees or limiting the number of trades.

Working with leverage
~~~~~~~~~~~~~~~~~~~~~

By default all positions traded will be long positions. You can allow
allocations in short positions by definining a nonzero limit on the total short
allcocations.  For example, to allow short selling up to 30% of the
portfolio value, you can do:

.. testcode:: mod

    import pandas as pd
    from gurobi_optimods.datasets import load_portfolio
    from gurobi_optimods.portfolio import MeanVariancePortfolio
    data = load_portfolio()
    Sigma = data.cov()
    mu = data.mean()
    gamma = 100.0
    mvp = MeanVariancePortfolio(Sigma, mu)
    x = mvp.efficient_portfolio(gamma, max_total_short=0.3)

.. testoutput:: mod
    :hide:

    ...
    Optimize a model with 42 rows, 50 columns and 110 nonzeros
    ...
    Model has 55 quadratic objective terms
    ...

With leverage allowed we now obtain an optimal portfolio with three short
positions, totaling to about 14% of the wealth:

.. doctest:: mod
    :options: +NORMALIZE_WHITESPACE

    >>> x
        AA    0.437482
        BB    0.020704
        CC   -0.080789
        DD    0.271877
        EE    0.019897
        FF   -0.029849
        GG    0.083466
        HH    0.240992
        II    0.066809
        JJ   -0.030588
    dtype: float64

    >>> x[x<0].sum()
    -0.14122620800822816

Restricting the number of trades
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is possible to compute an optimal portfolio under the additional restriction
that only a limited number of positions can be traded.  This can be set through
the ``max_trades`` keyword parameter.  For the example above, restricting the
total number of trades to three we get the following optimal portfolio.

.. testcode:: mod

    import pandas as pd

    from gurobi_optimods.datasets import load_portfolio
    from gurobi_optimods.portfolio import MeanVariancePortfolio

    data = load_portfolio()
    Sigma = data.cov()
    mu = data.mean()
    gamma = 100.0

    mvp = MeanVariancePortfolio(Sigma, mu)
    x = mvp.efficient_portfolio(gamma, max_trades=3)

.. testoutput:: mod
    :hide:

    ...
    Optimize a model with 43 rows, 50 columns and 120 nonzeros
    ...
    Model has 55 quadratic objective terms
    ...

The returned solution now suggests to trade only the assets AA, DD, HH.

.. doctest:: mod
    :options: +NORMALIZE_WHITESPACE

    >>> x
	AA    0.482084
	BB    0.000000
	CC    0.000000
	DD    0.282683
	EE    0.000000
	FF    0.000000
	GG    0.000000
	HH    0.235233
	II    0.000000
	JJ    0.000000
    dtype: float64

Transaction fees
~~~~~~~~~~~~~~~~

In order to define fixed costs per transaction suggested by the optimal
portfolio :math:`x`, you can use the ``fees_buy`` (for long positions) and
``fees_sell`` (for short positions) keyword parameters.

.. testcode:: mod

    import pandas as pd

    from gurobi_optimods.datasets import load_portfolio
    from gurobi_optimods.portfolio import MeanVariancePortfolio

    data = load_portfolio()
    Sigma = data.cov()
    mu = data.mean()
    gamma = 100.0

    mvp = MeanVariancePortfolio(Sigma, mu)
    x = mvp.efficient_portfolio(gamma, fees_buy=0.005)

.. testoutput:: mod
    :hide:

    ...
    Optimize a model with 42 rows, 50 columns and 110 nonzeros
    ...
    Model has 55 quadratic objective terms
    ...

Note that the ``fees_buy`` parameter designates the transaction cost
*relative* to the total portfolio value.  In the above example we used
the value 0.005, meaning that each trasaction has a fixed-cost of 0.5%
of the total portfolio value.

All transaction fees are assumed to be covered by the portfolio itself,
thus reducing the total sum of the returned optimal portfolio:

.. doctest:: mod
    :options: +NORMALIZE_WHITESPACE

    >>> x.sum()
    0.95

Minimum position constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A minimum fraction of investment can be enforced upon each individual position,
preventing trades at negligible volume.  Use the keyword parameters
``min_long`` and ``min_short`` to set a thresholds for trading long and short
positions.  For example, here we enforce that at least 5% of the wealth are
allocated to each trade:

.. testcode:: mod

    import pandas as pd
    from gurobi_optimods.datasets import load_portfolio
    from gurobi_optimods.portfolio import MeanVariancePortfolio
    data = load_portfolio()
    Sigma = data.cov()
    mu = data.mean()
    gamma = 100.0
    mvp = MeanVariancePortfolio(Sigma, mu)
    x_plain = mvp.efficient_portfolio(gamma, max_total_short=0.3)
    x_minpos = mvp.efficient_portfolio(gamma, max_total_short=0.3, min_long=0.05, min_short=0.05)

.. testoutput:: mod
    :hide:

    ...
    Optimize a model with 42 rows, 50 columns and 110 nonzeros
    ...
    Model has 55 quadratic objective terms
    ...
    Optimize a model with 62 rows, 50 columns and 150 nonzeros
    ...
    Model has 55 quadratic objective terms
    ...

Comparing the two portfolios ``x_plain``, which has no minimum position
constraints set with ``x_minpos``, which defines these constraints, we see that
the latter portfolio is free of "tiny" transactions.

.. doctest:: mod
    :options: +NORMALIZE_WHITESPACE

    >>> pd.concat([x_plain, x_minpos], keys=["plain", "minpos"], axis=1)
           plain    minpos
    AA  0.437482  0.431366
    BB  0.020704  0.000000
    CC -0.080789 -0.070755
    DD  0.271877  0.284046
    EE  0.019897  0.000000
    FF -0.029849 -0.050000
    GG  0.083466  0.097149
    HH  0.240992  0.244677
    II  0.066809  0.063517
    JJ -0.030588  0.000000
