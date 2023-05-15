QUBO: Quadratic Unconstrained Binary Optimization
=================================================

A quadratic unconstrained binary optimization (QUBO) problem is a combinatorial
optimization problem with applications from finance and economics
to machine learning (:footcite:t:`kochenberger2014ubqp`).
A QUBO model is a special case of a mixed-integer quadratic program (MIQP) that has
only binary decision variables and no constraints.

QUBO is :math:`\mathcal{NP}`-hard and for many classical problems such as
maximum cut, graph coloring, and the partition problem,
embeddings into QUBO have been formulated (:footcite:t:`glover2018tutorial`).

QUBO models are also related to Ising models and used in adiabatic
quantum computation, where a physical process called quantum annealing solves the models.
As a consequence, there is a trend to convert arbitrary combinatorial
optimization problems to QUBO models (which is a non-trivial procedure) to solve them
finally with a quantum annealer.
However, in many cases it is more efficient to directly solve the original problem
with classical optimization methods and avoid the translation to QUBO.


Problem Specification
---------------------

.. tabs::

    .. tab:: Domain-Specific Description

        We are given a set :math:`I` of :math:`n` items,
        weights :math:`q_i \in \mathbb{R}`
        for each single item :math:`i \in I`,
        and weights :math:`q_{ij} \in \mathbb{R}` for each pair of distinct items
        :math:`i,j \in I,~ i \neq j`.

        The objective is to find a subset :math:`S^* \subseteq I` of items such that the sum of weights
        of all single items in :math:`S^*` and all pairs of distinct items in :math:`S^*` is
        minimal, i.e.,

        .. math::
            S^* = \arg \min_{S \subseteq I} \sum_{i \in S} q_i + \sum_{i,j \in S,~ i \neq j} q_{ij}

        We arrange the weights in an upper triangular matrix :math:`Q \in \mathbb{R}^{n \times n}`
        where entry :math:`q_{ij}` with :math:`i < j` is the weight for item pair :math:`i,j`, and
        entry :math:`q_{ii} = q_i` is the weight for single item :math:`i`.

        Note that the input matrix does not necessarily need to be in upper triangular format.
        We accept matrices :math:`Q'` that are populated in an arbitrary way and accumulate symmetric entries,
        i.e., :math:`q_{ij} = q'_{ij} + q'_{ji}` for all item pairs :math:`i,j` with :math:`i < j`.

    .. tab:: Optimization Model

        We define a binary decision vector :math:`x \in \{0,1\}^n` with variables

        .. math::
            x_i = \begin{cases}
                1 & \text{if item}\,i \in S^*\\
                0 & \text{otherwise.} \\
            \end{cases}

        The QUBO model is then given as:

        .. math::
            \begin{align}
            \max \quad        & x' Q x = \sum_{i \in I} \sum_{j \in I} q_{ij} x_i x_j & \\
            \mbox{s.t.} \quad & x \in \{0, 1\}^n &
            \end{align}

        Note that weights :math:`q_i = q_{ii}` for single items :math:`i \in I` are
        correctly considered in the objective function since :math:`x_i x_i = x_i` holds
        for binary variables.

The input data consisting of the item (pair) weights is defined as a matrix (see the
domain-specific description), either as a NumPy array :class:`numpy.ndarray`
or as a SciPy sparse matrix :class:`scipy.sparse`.

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

    result = solve_qubo(Q, output=True, log_file="gurobi.log")

.. testoutput:: qubo
    :hide:

    ...
    New QUBO solution found with objective -4.0

The model is solved as an MIQP by Gurobi.

.. collapse:: View Gurobi Logs

    .. code-block:: text

        Gurobi 10.0.1 (linux64) logging started Fri Apr 28 17:24:54 2023

        Set parameter LogFile to value "gurobi.log"
        Gurobi Optimizer version 10.0.1 build v10.0.1rc0 (linux64)

        CPU model: Intel(R) Core(TM) i7-8665U CPU @ 1.90GHz, instruction set [SSE2|AVX|AVX2]
        Thread count: 4 physical cores, 8 logical processors, using up to 8 threads

        Optimize a model with 0 rows, 3 columns and 0 nonzeros
        Model fingerprint: 0x0d77f9fa
        Model has 5 quadratic objective terms
        Variable types: 0 continuous, 3 integer (3 binary)
        Coefficient statistics:
        Matrix range     [0e+00, 0e+00]
        Objective range  [0e+00, 0e+00]
        QObjective range [2e+00, 6e+00]
        Bounds range     [1e+00, 1e+00]
        RHS range        [0e+00, 0e+00]
        Found heuristic solution: objective 0.0000000
        Found heuristic solution: objective -1.0000000
        Found heuristic solution: objective -4.0000000
        Presolve removed 0 rows and 3 columns
        Presolve time: 0.00s
        Presolve: All rows and columns removed

        Explored 0 nodes (0 simplex iterations) in 0.00 seconds (0.00 work units)
        Thread count was 1 (of 8 available processors)

        Solution count 3: -4 -1 0
        No other solutions better than -4

        Optimal solution found (tolerance 1.00e-04)
        Best objective -4.000000000000e+00, best bound -4.000000000000e+00, gap 0.0000%

        User-callback calls 84, time in user-callback 0.00 sec


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
