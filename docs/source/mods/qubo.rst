Quadratic Unconstrained Binary Optimization (QUBO)
==================================================

Many :math:`\mathcal{NP}`-hard discrete optimization problems that naturally
arise in application fields such as finance, energy, healthcare, and machine learning,
can be mapped to a quadratically unconstrained binary optimization (QUBO) problem
(:footcite:t:`kochenberger2014ubqp`), or equivalently to an Ising Hamiltonian
(:footcite:t:`lucas2014ising`).
Examples of such problems are max-cut, graph colouring, partitioning, and maximum
independent set. A QUBO model is a mixed-integer quadratic program (MIQP) where
the decision variables are restricted to take binary values and there are no
explicit constraints.

Since solving optimization problems typically requires significant computing resources,
research in the development of novel computing architectures designed to solve QUBO problems
has flourished in recent years (:footcite:t:`johnson2011quantum,matsubara2018ising,farhi2014quantum`).
Such devices are fast, general-purpose, and power efficient. However, having to
transform an optimization problem into a QUBO problem which is not as expressive and rich
as MIP has significant drawbacks. It can make the problem significantly harder by
penalizing the constraints into the objective function and expanding the solution space
including both feasible and infeasible solutions. This remains to be an obstacle to the
practical usage of special-purpose hardware in solving practically relevant optimization
problems. Therefore, it is more efficient to avoid QUBO translations, if possible, and
solve the problem as a MIP where the constraints are directly represented.


Problem Specification
---------------------

We are given a set :math:`I` of :math:`n` items, weights :math:`q_i \in
\mathbb{R}` for each single item :math:`i \in I`, and weights :math:`q_{ij} \in
\mathbb{R}` for each pair of distinct items :math:`i,j \in I,~ i \neq j`.

The objective is to find a subset :math:`S^* \subseteq I` of items such that the
sum of weights of all single items in :math:`S^*` and all pairs of distinct
items in :math:`S^*` is minimal, i.e.,

.. math::
    S^* = \arg \min_{S \subseteq I} \sum_{i \in S} q_i + \sum_{i,j \in S,~ i \neq j} q_{ij}

We arrange the weights in an upper triangular matrix :math:`Q \in \mathbb{R}^{n
\times n}` where entry :math:`q_{ij}` with :math:`i < j` is the weight for item
pair :math:`i,j`, and entry :math:`q_{ii} = q_i` is the weight for single item
:math:`i`.

Note that the input matrix does not necessarily need to be in upper triangular
format. We accept matrices :math:`Q'` that are populated in an arbitrary way and
accumulate symmetric entries, i.e., :math:`q_{ij} = q'_{ij} + q'_{ji}` for all
item pairs :math:`i,j` with :math:`i < j`.

.. dropdown:: Background: Optimization Model

    This Mod is implemented by formulating the QUBO problem as a Binary
    Quadratic Program (BQP). To do so, we define a binary decision vector
    :math:`x \in \{0,1\}^n` with variables

    .. math::
        x_i = \begin{cases}
            1 & \text{if item}\,i \in S^*\\
            0 & \text{otherwise.} \\
        \end{cases}

    The BQP is then formulated as:

    .. math::
        \begin{align}
        \max \quad        & x' Q x = \sum_{i \in I} \sum_{j \in I} q_{ij} x_i x_j & \\
        \mbox{s.t.} \quad & x \in \{0, 1\}^n &
        \end{align}

    Note that weights :math:`q_i = q_{ii}` for single items :math:`i \in I` are
    correctly considered in the objective function since :math:`x_i x_i = x_i` holds
    for binary variables.

The input data consisting of the item (pair) weights is defined as a matrix (see the
description), either as a NumPy array :class:`~numpy.ndarray`
or as a SciPy sparse matrix :class:`~scipy.sparse.spmatrix`.

Code
----

The example below solves a QUBO problem instance based on 3 items
with single-item weights :math:`q_1 = 0,~ q_2 = -3,~ q_3 = 2`, and
item-pair weights :math:`q_{12} = -1,~ q_{13} = -2,~ q_{23} = 3`,
resulting in the following item weight matrix:

.. math::
    Q = \begin{pmatrix}
    0 & -1 & -2\\
    0 & -3 & 3\\
    0 & 0 & 2
    \end{pmatrix}

We use a NumPy array to represent matrix :math:`Q` (and alternatively we show the
definition as a SciPy sparse matrix in a comment).

.. testcode:: qubo

    import numpy as np
    import scipy.sparse as sp
    from gurobi_optimods.qubo import solve_qubo

    Q = np.array([[0, -1, -2], [0, -3, 3], [0, 0, 2]])

    # weights = [-3, 2, -1, -2, 3]
    # row = [1, 2, 0, 0, 1]
    # col = [1, 2, 1, 2, 2]
    # Q = sp.coo_matrix((weights, (row, col)), shape=(3, 3))

    result = solve_qubo(Q)

.. testoutput:: qubo
    :hide:

    ...
    New QUBO solution found with objective -4.0

Solution
--------

The returned result is a data class containing the objective value and
the solution itself as a NumPy ndarray.

.. doctest:: qubo
    :options: +NORMALIZE_WHITESPACE

    >>> result
    QuboResult(solution=array([1., 1., 0.]), objective_value=-4.0)
    >>> result.objective_value
    -4.0
    >>> result.solution
    array([1., 1., 0.])

.. footbibliography::
