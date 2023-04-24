QUBO: Quadratic Unconstrained Binary Optimization
=================================================

A quadratic unconstrained binary optimization (QUBO) problem is a combinatorial
optimization problem with a wide range of applications from finance and economics
to machine learning (:footcite:t:`kochenberger2014ubqp`).
A QUBO model is a special case of a mixed-integer quadratic program (MIQP) that has
only binary decision variables and no constraints.

QUBO is :math:`\mathcal{NP}`-hard and for many classical problems such as
maximum cut, graph coloring, and the partition problem,
embeddings into QUBO have been formulated (:footcite:t:`glover2018tutorial`).

QUBO models are also related to Ising models and therefore used in adiabatic
quantum computation, where a physical process called quantum annealing solves them.
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

        The objective is to find a subset :math:`S \subseteq I` of items such that the sum of weights
        of all single items in :math:`S` and all pairs of distinct items in :math:`S` is
        minimal, i.e.,

        .. math::
            S = \arg \min_{S' \subseteq I} \sum_{i \in S'} q_i + \sum_{i,j \in S',~ i \neq j} q_{ij}

        We arrange the weights in an upper triangular matrix :math:`Q \in \mathbb{R}^{n \times n}`
        where entry :math:`q_{ij}` with :math:`i < j` is the weight for item pair :math:`i,j`, and
        entry :math:`q_{ii} = q_i` is the weight for single item :math:`i`.

        Note that the input matrix does not necessarily need to be in upper triangular format.
        We accept matrices :math:`Q'` that are populated in an arbitrary way and accumulate symmetric entries,
        i.e., :math:`q_{ij} = q'_{ij} + q'_{ji}` for all items :math:`i < j`.

    .. tab:: Optimization Model

        We define a binary decision vector :math:`x \in \{0,1\}^n` with variables

        .. math::
            x_i = \begin{cases}
                1 & \text{if item}\,i \in S\\
                0 & \text{otherwise.} \\
            \end{cases}

        The QUBO model is then given as:

        .. math::
            \begin{align}
            \max \quad        & x' Q x = \sum_{i \in I} \sum_{j \in I} q_{ij} x_i x_j & \\
            \mbox{s.t.} \quad & x \in \{0, 1\}^n &
            \end{align}


.. footbibliography::
