Cardinality Constrained Regression
==================================

``NOT IMPLEMENTED YET``

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

    .. tab:: Mathematical Model

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

.. literalinclude:: ../../../examples/card_regression.py
