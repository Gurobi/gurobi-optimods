ACOPF
=====


Problem Specification
---------------------

.. tabs::

    .. tab:: Domain-Specific Description

        In the standard ACOPF problem we set a minimum-cost operating point (complex voltages at all buses, active and reactive power generation at each generator) for a power system so that the correct amount of active and reactive power is delivered to each bus and all branch limits are satisfied; both using the AC power flow laws.  Under this model, cost is incurred at the generators.  At each generator we are given a convex quadratic or piecewise-linear (convex) cost function associated with active power generation.

    .. tab:: Optimization Model
       The following variables are used in the ACOPF problem.
       1. For each bus :math:`k` we have a complex voltage :math:`V_k = |V_k| e^{j \theta_k}` (where :math:`j = \sqrt{-1}`).
       2. For each branch :math:`km` we have the two complex variables :math:`S_{km} \, \text{and} \, S_{mk}`, the complex power injected into the branch at :math:`k \, \text{and at} \, m`, respectively. (Recall that branches are given as ordered pairs).

       3. For each generator :math:`i` there is a (complex) generation :math:`P^{g}_i + j Q^{g}_i`.

       There are five classes of constraints: branch AC power definitions, flow balance, branch limits, generator limits and voltage bounds.

       The first set of constraints embody AC power flow laws. For each branch :math:`km` we have

       .. math::

           S _{km} =  V_k \left[ Y \left( \begin{array}{r}
	   V_{k} \\
	   V_{m}
	   \end{array} \right) \right]^*_k \quad \text{and} \quad S _{mk} =  V_m \left[ Y \left( \begin{array}{r}
	   V_{k} \\
	   V_{m}
	   \end{array} \right) \right]^*_m,

       where * indicates complex conjugate.  These are the only nonconvex expressions in the entire model.  Elsewhere we describe alternative representations of these constraints using different notation.

       The next set of constraints represent flow balance.  For each branch :math:`k` we have

       .. math::

	  \sum_{i \in G(k)} (P^{g}_i + jQ^{g}_i) \ - \ (P^{d}_k + j Q^{d}_k) \quad = \quad \sum_{km \in \delta^+(k)} S_{km} \ + \ \sum_{mk \in \delta^-(k)}S_{km}.


       In this expression, the left-hand side is the net excess of generation over load at bus :math:`k`.  The right-hand side is the total power injected into the grid at bus :math:`k`.

       Next we have branch limits.  For each branch :math:`km`,

       .. math::

	  \max\{ |S_{km}|, |S_{mk}| \} \ \le \, L_{km}.

       Note: there are alternative versions of this particular constraint.  In a power system under normal (non-stressed) operations this constraint is slack for all but a small number of branches.

       Finally we have generator limits and voltage bounds.

       .. math::

	  P^{\min}_i \, \le \, P^{g}_i \, \le \, P^{\max}_i \quad \text{and} \quad  Q^{\min}_i \, \le \, Q^{g}_i \, \le \, Q^{\max}_i, \ \text{for each generator} \, i,

       and

       .. math::
	  L_k \, \le \, |V_k| \, \le \, U_k, \ \text{for each bus} \, k.


       The objective to optimize is the following cost expression:


       .. math::
	  \min \sum_{\text{generators} \, i} F_i(P^g_i).


       Above we provided an expression defining the complex power flow :math:`S_{km}`in terms of the admittance matrix math::`Y`.  Denoting

       .. math::
	  \theta_{km} = \theta_k - \theta_m - \phi_{km},

       The real and imaginary components of :math:`S_{km}`

       .. math::
	  P_{km} \ = \ G_{kk} |V_k|^2 + G_{km}|V_k||V_m| \cos(\theta_{km}) + B_{km}|V_k||V_m| \sin(\theta_{km})

       and

       .. math::
	  Q_{km} \ = \ -B_{kk} |V_k|^2 - B_{km}|V_k||V_m| \cos(\theta_{km}) + G_{km}|V_k||V_m| \sin(\theta_{km})


       Likewise, in terms of :math:`S_{mk}` we have

       .. math::
	  P_{mk} \ = \ G_{mm} |V_m|^2 + G_{mk}|V_k||V_m| \cos(\theta_{km}) - B_{mk}|V_k||V_m| \sin(\theta_{km})

       and

       .. math::
	  Q_{mk} \ = \ -B_{mm} |V_m|^2 - B_{mk}|V_k||V_m| \cos(\theta_{km}) - G_{mk}|V_k||V_m| \sin(\theta_{km})



       By introducing the auxiliary variables

       .. math::
	  v_{k}^{(2)} \ = \ |V_k|^2 \ \text{for every bus} \, k

       and

       .. math::
	  c_{km} \ = \ |V_k||V_m| \cos(\theta_{km}), \quad s_{km} \ = \ |V_k||V_m| \sin(\theta_{km}) \ \text{for every branch} \, km

       the power flow quantities can be rewritten as

       .. math::
	  P_{km} \ = \ G_{kk} v_k^{(2)} + G_{km} c_{km} + B_{km} s_{km}

       and

       .. math::
	  Q_{km} \ = \ -B_{kk} v_k^{(2)}  - B_{km} c_{km} + G_{km} s_{km}.

       respectively.

    .. tab:: Jabr relaxation

       We can obtain a SOC (second-order cone) relaxation of the ACOPF formulation
       by introducing the :math:`v^{(2)}_k, \ c_{km} \ \text{and} \ s_{km}` auxiliary variables, removing the nonconvex definitions of such variables (which involve cosines and sines) and adding the rotated cone constraints

       .. math::
	  c_{km}^2 \ + \ s_{km}^2 \ \le \ v_k^{(2)} v_m^{(2)} \ \text{for every branch} \, km.


       The resulting relaxation can prove very tight, though challenging in large cases.

    .. tab:: QCQP

       ACOPF can be reformulated as a QCQP (quadratically-constrained quadratic program) by performing two reformulation steps.

       1. For each bus :math:`k`, introduce the real variables :math:`e_k \, \text{and} \, f_k` and set :math:`v^{(2)}_k \, = \, e_k^2 + f_k^2`.

       2. For each branch :math:`km`, set :math:`c_{km} = e_k e_m + f_k f_m` and :math:`s_{km} = -e_k f_m + f_k e_m`.

       These constraints render an exact reformulation rendering the problem as a QCQP, i.e., one that removes the sines and cosines from the formulation.  We remind the reader that, for example, we are writing the active power flow injected at bus :math:`k` on branch :math:`km` through the constraints

       .. math::
	  P_{km} \ = \ G_{kk} v_k^{(2)} + G_{km} c_{km} + B_{km} s_{km},

       (and similarly with :math:`Q_{km}, \, P_{mk}, \,\text{and} \,Q_{mk}`).








Code
----

Self contained code example to run the mod from an example dataset. Example
datasets should bd included in the ``gurobi_optimods.datasets`` module for
easy access by users.

.. testcode:: mod

    from gurobi_optimods.opf import solve_opf_model, read_case_from_mat_file
    from gurobi_optimods.datasets import load_caseopfmat

    settings = {"doac": True, "use_ef": True}
    # load path to case file
    casefile = load_caseopfmat("9")
    # read case file and return a case dictionary
    case = read_case_from_mat_file(casefile)
    # solve opf model and return a solution and the final objective value
    solution = solve_opf_model(settings, case)

..  A snippet of the Gurobi log output here won't show in the rendered page,
    but serves as a doctest to make sure the code example runs. The ... lines
    are meaningful here, they will match anything in the output test.

.. testoutput:: mod
    :hide:

    ...
    Optimize a model with 73 rows, 107 columns and 208 nonzeros
    ...
    Optimal solution found (tolerance 1.00e-03)
    ...

The model is solved as an LP/MIP/QP by Gurobi.

..  You can include the full Gurobi log output here for the curious reader.
    It will be visible as a collapsible section.

.. collapse:: View Gurobi Logs

    .. code-block:: text

        Gurobi Optimizer version 9.5.1 build v9.5.1rc2 (mac64[x86])
        Optimize a model with ...
        Best obj ... Best bound ...

|

Solution
--------

Show the solution. One way is to use doctests to display simple shell outputs
(see the workforce example). This can be done simply by pasting outputs
directly from a python shell. Another option is to include and display figures
(see the graph matching examples).

.. doctest:: mod
    :options: +NORMALIZE_WHITESPACE

    >>>
