Cardinality Constrained Regression
===================================

Cardinality Constrained Regression (CCR) provides a powerful approach to feature
selection by explicitly limiting the number of non-zero coefficients in a linear
model. Unlike regularization methods such as Lasso, which apply soft penalties to
encourage sparsity, CCR uses mixed-integer optimization to guarantee that exactly
:math:`k` features are selected. This makes CCR particularly valuable in scenarios
where hard constraints on model complexity are required for interpretability,
operational feasibility, or cost management.

Motivation: When Hard Constraints Matter
-----------------------------------------

Consider a hospital developing a risk scoring system to predict patient readmission
within 30 days of discharge. The challenge extends beyond predictive accuracy:

* **Operational constraints**: Each clinical test costs money and time. Lab work,
  imaging studies, and specialist consultations must be budgeted. The hospital can
  only afford to systematically collect :math:`k=5` predictive markers for every
  discharged patient.

* **Clinical interpretability**: Healthcare providers need to understand and trust
  the model. A scoring system with exactly 5 predictors can be memorized, explained
  to patients, and applied at bedside without computer assistance.

* **Regulatory compliance**: Many healthcare systems require transparent, auditable
  decision-making processes. A model with a guaranteed small number of features is
  easier to document and defend.

In such scenarios, traditional regularization methods fall short. Lasso might
select 15-20 features with varying coefficient magnitudes, and the exact number
depends on the tuning parameter :math:`\alpha`. By contrast, CCR provides a hard
guarantee: exactly :math:`k` features will be non-zero, allowing stakeholders to
make firm operational and budgetary commitments.

Similar applications arise in financial credit scoring for underbanked populations,
environmental monitoring with sensor budget constraints, and manufacturing quality
control with limited inspection resources. In each case, the combinatorial nature
of feature selection—choosing the optimal subset of exactly :math:`k` from :math:`n`
features—is naturally expressed as a mixed-integer optimization problem.

Problem Specification
---------------------

``CardinalityConstrainedRegression`` fits a linear model with coefficients
:math:`w` to minimize the sum of squared errors, subject to a hard constraint
on the number of non-zero coefficients. This is formulated as a mixed-integer
quadratic program (MIQP):

.. math::

    \begin{alignat}{2}
    \min_w \quad        & \lVert Xw - y \rVert_2^2 \\
    \mbox{s.t.} \quad   & \lVert w \rVert_0 \le k
    \end{alignat}

Here, :math:`\lVert w \rVert_0` denotes the number of non-zero entries in
:math:`w` (the L0 "norm"), and :math:`k` is the cardinality constraint
specifying the maximum number of features to include in the model.

.. dropdown:: Background: Mathematical Formulation

    The optimization problem is solved by Gurobi as a mixed-integer quadratic
    program. Let :math:`I` denote the set of observations and :math:`J` the set
    of features. Response values :math:`y_i` are predicted from predictor values
    :math:`x_{ij}` by fitting coefficients :math:`w_j`. Binary variables
    :math:`b_j` track whether each coefficient is non-zero, and the cardinality
    constraint limits the sum of these indicators.

    .. math::

        \begin{alignat}{2}
        \min \quad        & \sum_i e_i^2 \\
        \mbox{s.t.} \quad & \sum_j w_j x_{ij} + e_i = y_i  \quad & \forall i \in I \\
                          & w_j \in \{0\} \cup [-M, M]     \quad & \forall j \in J \\
                          & w_j = 0 \vee b_j = 1           \quad & \forall j \in J \\
                          & \sum_j (1 - b_j) \le k \\
                          & e_i \,\, \text{free}           \quad & \forall i \in I \\
                          & w_j \,\, \text{free}           \quad & \forall j \in J \\
                          & b_j \in \{0, 1\}               \quad & \forall j \in J
        \end{alignat}

    The relationship :math:`w_j = 0 \vee b_j = 1` is enforced using Special
    Ordered Sets (SOS) constraints, which Gurobi handles efficiently.

Example: Hospital Readmission Risk Scoring
-------------------------------------------

This example demonstrates how cardinality constrained regression can be used to
develop an interpretable clinical risk score. We use the diabetes dataset from
scikit-learn as a proxy for a healthcare prediction task. The goal is to build
a model that uses exactly :math:`k=5` clinical features to predict patient
outcomes, ensuring the model remains simple enough for clinical deployment.

Basic Usage
^^^^^^^^^^^

The API follows scikit-learn conventions with ``fit`` and ``predict`` methods:

.. testcode:: card_regression

    from sklearn import datasets
    from sklearn.model_selection import train_test_split

    from gurobi_optimods.regression import CardinalityConstrainedRegression

    # Load the diabetes dataset
    diabetes = datasets.load_diabetes()

    # Split data for fit assessment
    X_train, X_test, y_train, y_test = train_test_split(
        diabetes["data"], diabetes["target"], random_state=42, test_size=0.2
    )

    # Create and fit model with exactly 5 non-zero coefficients
    ccr = CardinalityConstrainedRegression(k=5)
    ccr.fit(X_train, y_train)
    y_pred = ccr.predict(X_test)

.. testoutput:: card_regression
    :hide:

    ...
    Optimize a model with ...
    ...
    Optimal solution found (tolerance 1.00e-04)
    Best objective ...

The model guarantees that exactly 5 features will have non-zero coefficients,
making it simple to implement in a clinical setting.

Comparison with Lasso Regularization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To understand the value of cardinality constraints, we compare CCR with Lasso
regression. While Lasso encourages sparsity through an L1 penalty, it does not
provide hard guarantees on the number of features selected.

.. testcode:: card_regression

    import pandas as pd
    import matplotlib.pyplot as plt
    from sklearn.linear_model import Lasso

    # Fit Lasso with a regularization parameter
    lasso = Lasso(alpha=1.0, random_state=42)
    lasso.fit(X_train, y_train)

    # Compare coefficients
    coefficients = pd.DataFrame(
        data={
            "CCR (k=5)": ccr.coef_,
            "Lasso (α=1.0)": lasso.coef_
        },
        index=diabetes["feature_names"]
    )

    # Visualize the selected features
    fig, ax = plt.subplots(figsize=(10, 5))
    coefficients.plot.bar(ax=ax)
    ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
    ax.set_ylabel("Coefficient Value")
    ax.set_xlabel("Clinical Feature")
    ax.set_title("Feature Selection: CCR vs Lasso")
    plt.tight_layout()

    # Count non-zero coefficients
    n_nonzero_ccr = (coefficients["CCR (k=5)"].abs() > 1e-6).sum()
    n_nonzero_lasso = (coefficients["Lasso (α=1.0)"].abs() > 1e-6).sum()
    print(f"CCR selected features: {n_nonzero_ccr}")
    print(f"Lasso selected features: {n_nonzero_lasso}")

.. testoutput:: card_regression
    :hide:

    ...
    CCR selected features: 5
    Lasso selected features: ...

The output demonstrates that CCR provides an exact guarantee on model complexity,
while Lasso's sparsity depends on the regularization parameter :math:`\alpha` and
may not achieve the desired level of simplicity.

Evaluating Model Performance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We assess both models' predictive accuracy and interpretability:

.. testcode:: card_regression

    from sklearn.metrics import mean_squared_error, r2_score

    # Predictions on test set
    y_pred_ccr = ccr.predict(X_test)
    y_pred_lasso = lasso.predict(X_test)

    # Compute metrics
    mse_ccr = mean_squared_error(y_test, y_pred_ccr)
    mse_lasso = mean_squared_error(y_test, y_pred_lasso)
    r2_ccr = r2_score(y_test, y_pred_ccr)
    r2_lasso = r2_score(y_test, y_pred_lasso)

    print(f"CCR    - MSE: {mse_ccr:.2f}, R²: {r2_ccr:.3f}")
    print(f"Lasso  - MSE: {mse_lasso:.2f}, R²: {r2_lasso:.3f}")

    # Show selected features for clinical interpretation
    print(f"\nSelected features for CCR model (k=5):")
    selected = coefficients[coefficients["CCR (k=5)"].abs() > 1e-6]["CCR (k=5)"]
    for feature, coef in selected.items():
        print(f"  {feature:10s}: {coef:8.3f}")

.. testoutput:: card_regression
    :hide:

    CCR    - MSE: ...
    Lasso  - MSE: ...
    ...
    Selected features for CCR model (k=5):
    ...

In a clinical deployment, the CCR model's limited feature set makes it practical
to implement. Healthcare providers can focus on collecting and monitoring exactly
5 key indicators, making the risk scoring system operationally feasible and
transparent to patients and regulators.

Choosing the Cardinality Parameter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The choice of :math:`k` involves a trade-off between model accuracy and
interpretability. We can evaluate this trade-off by fitting models with different
cardinality constraints:

.. testcode:: card_regression

    import numpy as np

    # Evaluate models with different cardinality constraints
    k_values = range(1, 11)
    results = []

    for k in k_values:
        ccr_k = CardinalityConstrainedRegression(k=k)
        ccr_k.fit(X_train, y_train, verbose=False)
        y_pred_train = ccr_k.predict(X_train)
        y_pred_test = ccr_k.predict(X_test)

        results.append({
            "k": k,
            "mse_train": mean_squared_error(y_train, y_pred_train),
            "mse_test": mean_squared_error(y_test, y_pred_test),
            "r2_test": r2_score(y_test, y_pred_test)
        })

    results_df = pd.DataFrame(results)

    # Visualize the trade-off
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

    ax1.plot(results_df["k"], results_df["mse_train"], 'o-', label="Train MSE")
    ax1.plot(results_df["k"], results_df["mse_test"], 's-', label="Test MSE")
    ax1.set_xlabel("Number of Features (k)")
    ax1.set_ylabel("Mean Squared Error")
    ax1.set_title("Model Complexity vs Error")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    ax2.plot(results_df["k"], results_df["r2_test"], 'o-', color='green')
    ax2.set_xlabel("Number of Features (k)")
    ax2.set_ylabel("R² Score (Test Set)")
    ax2.set_title("Predictive Performance vs Complexity")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    # Find the "elbow" point
    best_k = results_df.loc[results_df["mse_test"].idxmin(), "k"]
    print(f"\nBest k for test set performance: {best_k}")

.. testoutput:: card_regression
    :hide:

    ...
    Best k for test set performance: ...

This analysis helps identify the optimal cardinality constraint that balances
model simplicity with predictive accuracy. In practice, domain experts (e.g.,
clinicians) can use this information alongside operational constraints to select
an appropriate value of :math:`k`.

When to Use Cardinality Constrained Regression
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CCR is particularly valuable when:

* **Hard operational constraints exist**: Budget limitations, data collection
  costs, or processing time require an exact number of features.

* **Interpretability is critical**: Stakeholders need to understand and trust
  the model's decisions (healthcare, lending, legal applications).

* **Regulatory requirements demand transparency**: Auditable models with
  clear decision logic are required.

* **Features have acquisition costs**: Each additional measurement incurs
  financial or time costs (lab tests, sensor deployments, manual inspections).

In contrast, Lasso may be preferable when approximate sparsity is sufficient and
the computational cost of mixed-integer optimization is prohibitive for very
high-dimensional problems (thousands of features).

Technical Details
-----------------

The model is solved as a mixed-integer quadratic program (MIQP) by Gurobi.
The binary variables introduced for feature selection are handled using
Special Ordered Sets (SOS) constraints, which Gurobi's presolve and branch-and-bound
algorithms exploit efficiently.

.. collapse:: View Gurobi logs (example output)

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
