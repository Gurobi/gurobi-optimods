Least Absolute Deviations Regression
====================================

Minimum sum of absolute errors (L1) regression is generally more robust than ordinary least squares (OLS, L2) in that it is more resistant to outliers in the response variable. The loss function can be expressed using linear program (LP), so fitting model coefficients is ideally suited to an LP solver.

The interface of this mod matches that of :code:`sklearn.linear_model.LinearRegression`. This example compares the coefficients found using L1 and L2 regression on the diabetes dataset.

- Comparison to sklearn OLS: L1 norm is not implemented in sklearn. There is in general no analytic solution and gradient decent is not effective? So we need an LP solver.
- More robust (i.e. less sensitive to outliers) than L2 norm regression (OLS). Because the error metric is linear, an increase in the deviation of an individual point has a less extreme effect.

Problem Specification
---------------------

Scikit-learn's documentation gives a general explanation of `Linear Models <https://scikit-learn.org/stable/modules/linear_model.html>`_. The distinction between this mod and the Ordinary Least Squares model from scikit-learn is the loss function.

.. tabs::

    .. tab:: Loss Function

        :code:`LADRegression` fits a linear model with coefficients :math:`w` to minimize the sum of absolute errors.

        .. math::

            \min_w \lvert Xw - y \rvert

    .. tab:: Optimization Model

        To model the L1 regression loss function using linear programming, we need to introduce a number of auxiliary variables. Here :math:`I` is the set of data points and :math:`J` the set of fields. Response values :math:`y_i` are predicted from predictor values :math:`x_{ij}` by fitting coefficients :math:`w_j`. To handle the absolute value, non-negative variables :math:`u_i` and :math:`v_i` are introduced.

        .. math::

            \begin{alignat}{2}
            \min \quad        & \sum_i u_i + v_i \\
            \mbox{s.t.} \quad & \sum_j w_j x_{ij} + u_i - v_i = y_i \quad & \forall i \in I \\
                              & u_i, v_i \ge 0                      \quad & \forall i \in I \\
                              & w_j \,\, \text{free}                \quad & \forall j \in J \\
            \end{alignat}

*TODO: add data examples here (ndarrays)*

Code
----

This mod implements the fit-predict interface of scikit-learn. The example below reads in the diabetes dataset from scikit-learn, performs a train-test split, fits the L1 regression model to the training data, and creates predictions for the testing data.

.. testcode:: lad_regression

    from sklearn import datasets
    from sklearn.model_selection import train_test_split

    from gurobi_optimods.regression import LADRegression

    # Load the diabetes dataset
    diabetes_X, diabetes_y = datasets.load_diabetes(return_X_y=True)

    # Split data for fit assessment
    X_train, X_test, y_train, y_test = train_test_split(
        diabetes_X, diabetes_y, random_state=42
    )

    # Create and fit parameterised model
    reg = LADRegression()
    reg.fit(X_train, y_train)
    y_pred = reg.predict(X_test)

.. testoutput:: lad_regression
    :hide:

    ...
    Optimize a model with 331 rows, 673 columns and 4303 nonzeros
    ...
    Optimal objective  1.44...

The model is solved as a linear program by Gurobi. Logs provided for interested parties:

.. collapse:: View Gurobi logs

    .. code-block:: text

        Gurobi Optimizer version 10.0.1 build v10.0.1rc0 (mac64[x86])

        CPU model: Intel(R) Core(TM) i5-1038NG7 CPU @ 2.00GHz
        Thread count: 4 physical cores, 8 logical processors, using up to 8 threads

        Optimize a model with 331 rows, 673 columns and 4303 nonzeros
        Model fingerprint: 0xb54fc171
        Coefficient statistics:
          Matrix range     [6e-05, 1e+00]
          Objective range  [1e+00, 1e+00]
          Bounds range     [0e+00, 0e+00]
          RHS range        [2e+01, 3e+02]
        Presolve time: 0.00s
        Presolved: 331 rows, 673 columns, 4303 nonzeros

        Iteration    Objective       Primal Inf.    Dual Inf.      Time
               0      handle free variables                          0s
             354    1.4473274e+04   0.000000e+00   0.000000e+00      0s

        Solved in 354 iterations and 0.01 seconds (0.01 work units)
        Optimal objective  1.447327363e+04

|

Solution
--------

Here we extract the coefficients of the fitted model and compare them with the coefficients found using OLS. Not a super informative plot at this stage...

.. testcode:: lad_regression

    import pandas as pd
    from sklearn.linear_model import LinearRegression
    ols = LinearRegression()
    ols.fit(X_train, y_train)
    pd.DataFrame(data={"OLS": ols.coef_, "L1": reg.coef_}).plot.bar()

.. image:: figures/reg_coeffs.png
  :width: 500
  :alt: Weighted matching result

To gasps of shock and awe, the L1 regression produces a *smaller mean absolute error* on the training set than the OLS model, while the OLS model does better in terms of mean squared error.

.. doctest:: lad_regression

    >>> from sklearn.metrics import mean_absolute_error, mean_squared_error
    >>> round(mean_absolute_error(y_train, reg.predict(X_train)), 2)
    43.73
    >>> round(mean_absolute_error(y_train, ols.predict(X_train)), 2)
    44.05
    >>> round(mean_squared_error(y_train, reg.predict(X_train)), 1)
    2960.7
    >>> round(mean_squared_error(y_train, ols.predict(X_train)), 1)
    2907.3

Interesting related reading
---------------------------

- L1 regression is more commonly referred to as LAD (least absolute deviations) in the literature. I should probably change this terminology.
- `sklego <https://scikit-lego.netlify.app/linear-models.html#Least-Absolute-Deviation-Regression>`_ has an LAD implementation
- `Statsmodels <https://www.statsmodels.org/dev/regression.html>`_ has a quantile regression implementation (and the docs claim $q=0.5$ is equivalent to LAD)
- :footcite:t:`birkes2011alternative`
    - Chapter 4 intro quote: The method of least absolute deviations was introduced almost 50 years before the method of least squares, in 1757 by Roger Joseph Boscovich. He devised the method as a way to reconcile inconsistent measurements for the purpose of estimating the shape of the earth. After Pierre Simon Laplace adopted the method 30 years later, it saw occasional use, but it was soon overshadowed by the method of least squares. The popularity of least squares was at least partly due to the relative simplicity of its computations and to the supporting theory that was developed for it by Gauss and Laplace. Today, computation is not such a limitation and theoretical foundations have been laid for a variety of alternative methods, including the method of least absolute deviations (LAD).
    - Chapter 9 quote: The strength of LAD estimation is its robustness with respect to the distribution of the response variable (although not with respect to the explanatory variables).
- :footcite:t:`bloomfield1980least`
    - Idea predates least squares, but the computations are more complex
    - The development of linear programming made this problem manageable

.. footbibliography::
