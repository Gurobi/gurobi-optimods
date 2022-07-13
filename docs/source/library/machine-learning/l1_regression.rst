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

    .. tab:: Mathematical Model

        To model the L1 regression loss function using linear programming, we need to introduce a number of auxiliary variables. Here :math:`I` is the set of data points and :math:`J` the set of fields. Response values :math:`y_i` are predicted from predictor values :math:`x^j_i` by fitting coefficients :math:`w^j`. To handle the absolute value, non-negative variables :math:`u_i` and :math:`v_i` are introduced.

        .. math::

            \begin{alignat}{2}
            \min \quad        & \sum_i u_i + v_i \\
            \mbox{s.t.} \quad & \sum_j w^j x^j_i + u_i - v_i = y_i \quad & \forall i \in I \\
                              & u_i, v_i \ge 0                     \quad & \forall i \in I \\
                              & w^j \,\, \text{free}               \quad & \forall j \in J \\
            \end{alignat}

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

Output from the predictive model, just like in sklearn.

.. testcode:: l1_regression
    :hide:

    import sys
    sys.path.append("examples")
    from l1_regression_nupstup import y_pred, y_test
    sys.path.pop()

.. testoutput:: l1_regression
    :hide:

    Gurobi Optimizer version ...
    Optimal objective  4.372590220e+01

.. testcode:: l1_regression

    # Assess error
    from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
    print("Mean squared error: %.2f" % mean_squared_error(y_test, y_pred))
    print("Mean absolute error: %.2f" % mean_absolute_error(y_test, y_pred))
    print("Coefficient of determination: %.2f" % r2_score(y_test, y_pred))

.. testoutput:: l1_regression

    Mean squared error: 2969.58
    Mean absolute error: 41.92
    Coefficient of determination: 0.46
