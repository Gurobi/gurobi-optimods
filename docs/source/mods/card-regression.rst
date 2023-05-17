Cardinality Constrained Regression
==================================

- Present this in contrast to :code:`sklearn.linear_model.Lasso`.
- Cardinality constrained regression provides a stricter approach to model regularization.
- Instead of including a weighted objective term which tries to avoid unnecessary large coefficients in a model, cardinality constrained regression enforced a limit on the number of non-zero coefficients.
- Result is a model with limited complexity: explainable AI, woohoo!

Problem Specification
---------------------

See sklearn `Lasso <https://scikit-learn.org/stable/modules/linear_model.html#lasso>`_ for general explanation.

.. tabs::

    .. tab:: Loss Function

        :code:`CardConstrainedRegression` fits a linear model with coefficients :math:`w` to minimize the sum of absolute errors.

        .. math::

            \begin{alignat}{2}
            \min_w \quad        & \lvert Xw - y \rvert \\
            \mbox{s.t.} \quad   & {\lvert x \rvert}_0 \le k \\
            \end{alignat}

    .. tab:: Optimization Model

        To model the L1 regression loss function using linear programming, we need to introduce a number of auxiliary variables. Here :math:`I` is the set of data points and :math:`J` the set of fields. Response values :math:`y_i` are predicted from predictor values :math:`x_{ij}` by fitting coefficients :math:`w_j`. To handle the absolute value, non-negative variables :math:`u_i` and :math:`v_i` are introduced. Additionally, binary variables :math:`b_i` track the number of non-zero coefficients.

        .. math::

            \begin{alignat}{2}
            \min \quad        & \sum_i u_i + v_i \\
            \mbox{s.t.} \quad & \sum_j w_j x_{ij} + u_i - v_i = y_i \quad & \forall i \in I \\
                              & -M b_j \le w_j \le M b_j           \quad & \forall j \in J \\
                              & \sum_j b_j \le k \\
                              & u_i, v_i \ge 0                     \quad & \forall i \in I \\
                              & w_j \,\, \text{free}               \quad & \forall j \in J \\
            \end{alignat}

Code
----

.. testcode:: card_regression

    from sklearn import datasets
    from sklearn.model_selection import train_test_split

    from gurobi_optimods.regression import CardinalityConstrainedRegression

    # Load the diabetes dataset
    diabetes_X, diabetes_y = datasets.load_diabetes(return_X_y=True)

    # Split data for fit assessment
    X_train, X_test, y_train, y_test = train_test_split(
        diabetes_X, diabetes_y, random_state=42
    )

    # Create and fit parameterised model, including an intercept
    # but with at most two non-zero coefficients.
    reg = CardinalityConstrainedRegression(k=2)
    reg.fit(X_train, y_train)
    y_pred = reg.predict(X_test)

.. testoutput:: card_regression
    :hide:

    ...
    Optimize a model with 332 rows, 352 columns and 3982 nonzeros
    ...
    Optimal solution found (tolerance 1.00e-04)
    Best objective 5.67...


The model is solved as a MIP by Gurobi. Logs provided for interested parties:

.. collapse:: View Gurobi logs

    .. code-block:: text

        Gurobi Optimizer version 10.0.1 build v10.0.1rc0 (mac64[x86])

        CPU model: Intel(R) Core(TM) i5-1038NG7 CPU @ 2.00GHz
        Thread count: 4 physical cores, 8 logical processors, using up to 8 threads

        Optimize a model with 332 rows, 352 columns and 3982 nonzeros
        Model fingerprint: 0x3e2736be
        Model has 331 quadratic objective terms
        Model has 10 SOS constraints
        Variable types: 342 continuous, 10 integer (10 binary)
        Coefficient statistics:
          Matrix range     [6e-05, 1e+00]
          Objective range  [0e+00, 0e+00]
          QObjective range [2e+00, 2e+00]
          Bounds range     [1e+00, 1e+00]
          RHS range        [8e+00, 3e+02]
        Presolve time: 0.00s
        Presolved: 332 rows, 352 columns, 3982 nonzeros
        Presolved model has 10 SOS constraint(s)
        Presolved model has 331 quadratic objective terms
        Variable types: 342 continuous, 10 integer (10 binary)

        Root relaxation: objective 4.581241e+06, 686 iterations, 0.02 seconds (0.02 work units)

            Nodes    |    Current Node    |     Objective Bounds      |     Work
         Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

             0     0 4581241.43    0    8          - 4581241.43      -     -    0s
        H    0     0                    6370343.1269 4581241.43  28.1%     -    0s
        H    0     0                    6220498.1405 4581241.43  26.4%     -    0s
             0     0 4584952.94    0    8 6220498.14 4584952.94  26.3%     -    0s
             0     2 4584952.94    0    8 6220498.14 4584952.94  26.3%     -    0s
        *    5     6               2    6173920.5974 4819424.92  21.9%  10.2    0s
        *   12    10               3    6078925.9446 4878409.47  19.7%   9.0    0s
        H   34     8                    6028487.5605 5253130.81  12.9%   7.9    0s
        *   35     8               6    5673934.6855 5253130.81  7.42%   8.1    0s

        Explored 57 nodes (1028 simplex iterations) in 0.12 seconds (0.14 work units)
        Thread count was 8 (of 8 available processors)

        Solution count 6: 5.67393e+06 6.02849e+06 6.07893e+06 ... 6.37034e+06

        Optimal solution found (tolerance 1.00e-04)
        Best objective 5.673934685517e+06, best bound 5.673934685517e+06, gap 0.0000%
