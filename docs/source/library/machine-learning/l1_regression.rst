L1 Norm Regression
==================

See sklearn LinearRegression (L2 norm, no regularization), Ridge (L2 norm, L2 regularization), Lasso (L2 norm, L1 regularization). This model uses Gurobi to implement *L1* norm (also with L1 regularization).

Problem Specification
---------------------

Standard linear regression uses ordinary least squares ...

.. tabs::

    .. tab:: Data Specification

        We match the sklearn APIs: provide X_train and y_train as normal.

        Give the regression objective function here, like sklearn does.

    .. tab:: Mathematical Model

        Show the full LP implementation needed to get around all the absolute values.

Code
----

Alternate between the code required to run the model from the store vs how to implement directly in gurobipy. If you use nupstup, all the gurobi internals are handled for you; users interact with the 'solver' by passing dataframes to a given spec and receiving a dataframe as output. If you instead peek under the hood and use gurobipy, you have more options to extend the model with additional constraints and data.

.. tabs::
    .. tab:: nupstup function

        .. literalinclude:: ../../../examples/l1_regression_nupstup.py
            :linenos:

    .. tab:: gurobipy model

        .. literalinclude:: ../../../examples/l1_regression_gurobipy.py
            :linenos:


Both codes construct the same model and give the same result. The model is solved as a linear program by Gurobi.

.. collapse:: View Gurobi logs

    .. code-block:: text

        Gurobi Optimizer version 9.5.1 build v9.5.1rc2 (mac64[x86])
        Thread count: 4 physical cores, 8 logical processors, using up to 8 threads
        Optimize a model with 351 rows, 683 columns and 4343 nonzeros
        Model fingerprint: 0xaab99046
        Coefficient statistics:
        Matrix range     [2e-04, 1e+00]
        Objective range  [3e-03, 1e-02]
        Bounds range     [0e+00, 0e+00]
        RHS range        [2e+01, 3e+02]
        Presolve time: 0.01s
        Presolved: 351 rows, 683 columns, 4343 nonzeros

        Iteration    Objective       Primal Inf.    Dual Inf.      Time
            0      handle free variables                          0s
            359    5.9234728e+01   0.000000e+00   0.000000e+00      0s

        Solved in 359 iterations and 0.02 seconds (0.02 work units)
        Optimal objective  5.923472777e+01

|

Solution
--------

Solution is a an output from the predictive model, just like in sklearn.

.. testcode:: workforce
    :hide:

    import sys
    sys.path.append("examples")
    from l1_regression_nupstup import y_pred, y_test
    sys.path.pop()

.. testoutput:: workforce
    :hide:

    Gurobi Optimizer version ...
    Optimal objective  5.920257736e+01

.. testcode:: workforce

    # Assess error
    from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
    print("Mean squared error: %.2f" % mean_squared_error(y_test, y_pred))
    print("Mean absolute error: %.2f" % mean_absolute_error(y_test, y_pred))
    print("Coefficient of determination: %.2f" % r2_score(y_test, y_pred))

.. testoutput:: workforce

    Mean squared error: 2930.86
    Mean absolute error: 43.76
    Coefficient of determination: 0.47
