L1 Regression
=============

- Minimum sum of absolute errors (L1) regression performs is more robust than ordinary least squares (OLS, L2) in that it is more resistant to outliers in the response variable.
- Expressed as a linear program, so ideally suited to Gurobi.
- Present this in contrast to :code:`sklearn.linear_model.LinearRegression` (find a good comparison dataset).
- This implementation matches the sklearn APIs, can be used as a drop-in replacement.

Problem Specification
---------------------

See sklearn `Linear Models <https://scikit-learn.org/stable/modules/linear_model.html>`_ for general explanation.

.. tabs::

    .. tab:: Loss Function

        :code:`L1Regression` fits a linear model with coefficients :math:`w` to minimize the sum of absolute errors.

        .. math::

            \min_w \lvert Xw - y \rvert

    .. tab:: Optimization Model

        To model the L1 regression loss function using linear programming, we need to introduce a number of auxiliary variables. Here :math:`I` is the set of data points and :math:`J` the set of fields. Response values :math:`y_i` are predicted from predictor values :math:`x_{ij}` by fitting coefficients :math:`w_j`. To handle the absolute value, non-negative variables :math:`u_i` and :math:`v_i` are introduced.

        .. math::

            \begin{alignat}{2}
            \min \quad        & \sum_i u_i + v_i \\
            \mbox{s.t.} \quad & \sum_j w_j x_{ij} + u_i - v_i = y_i \quad & \forall i \in I \\
                              & u_i, v_i \ge 0                     \quad & \forall i \in I \\
                              & w_j \,\, \text{free}               \quad & \forall j \in J \\
            \end{alignat}

Code
----

Show the code required to run the model from the store. All the gurobi internals are handled for you; users interact with the 'solver' by passing dataframes to a given spec and receiving a dataframe as output.

.. literalinclude:: ../../../examples/l1_regression.py

The model is solved as a linear program by Gurobi.

.. collapse:: View Gurobi logs

    .. code-block:: text

        Gurobi Optimizer version 9.5.1 build v9.5.1rc2 (mac64[x86])
        Thread count: 4 physical cores, 8 logical processors, using up to 8 threads
        Optimize a model with 331 rows, 673 columns and 4303 nonzeros
        Model fingerprint: 0x6983ca17
        Coefficient statistics:
        Matrix range     [6e-05, 1e+00]
        Objective range  [3e-03, 3e-03]
        Bounds range     [0e+00, 0e+00]
        RHS range        [2e+01, 3e+02]
        Presolve time: 0.00s
        Presolved: 331 rows, 673 columns, 4303 nonzeros

        Iteration    Objective       Primal Inf.    Dual Inf.      Time
            0      handle free variables                          0s
            354    4.3725902e+01   0.000000e+00   0.000000e+00      0s

        Solved in 354 iterations and 0.01 seconds (0.01 work units)
        Optimal objective  4.372590220e+01

|

Solution
--------

.. testcode:: l1_regression
    :hide:

    from examples.l1_regression import reg, y_pred, y_test

.. testoutput:: l1_regression
    :hide:

    ...
    Optimal objective  4.372590220e+01

Properties of the predictive model, just like in sklearn.

.. doctest:: l1_regression

    >>> reg.coef_
    array([  16.7152629 , -306.19230544,  454.36833914,  508.02507763,
           -990.07434864,  414.38167986,  260.18885417,  483.00952994,
            678.56792495,   14.56067715])
    >>> reg.intercept_
    151.61357348161457

Output from the predictive model, just like in sklearn.

.. doctest:: l1_regression

    >>> from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
    >>> mean_squared_error(y_test, y_pred)
    2969.577566715166
    >>> mean_absolute_error(y_test, y_pred)
    41.9166462209382
    >>> r2_score(y_test, y_pred)
    0.4629757409105141
