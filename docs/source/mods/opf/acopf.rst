ACOPF
=====

Prior to specifying the mathematical structure of the problem it is important to present basic concepts and terminology used by the power systems community.

Basic Terminology
---------------------
 (to be added: buses, branches, generators, voltage, current, power)

Problem Specification
---------------------

We will first provide a definition of the standard problem.  In modern versions described below, network devices (FACTS, phase-shifters, impedance corrections are also incorporated into the model.

.. tabs::

    .. tab:: Domain-Specific Description

        In the standard ACOPF problem we set a minimum-cost operating point (complex voltages at all buses, active and reactive power generation at each generator) for a power system so that the correct amount of active and reactive power is delivered to each bus and all branch limits are satisfied; both using the AC power flow laws.  Under this model, cost is incurred at the generators.  At each generator we are given a convex quadratic or piecewise-linear (convex) cost function associated with active power generation.

    .. tab:: Optimization Model

As input data for the problem we have a power system ('power grid' for non-engineers) consisting of a network of buses and branches. Each branch is given as
an ordered pair :math:`km`, where :math:`k` is the "from" bus and :math:`m` is the "to" bus.

1. For each branch :math:`km` we have a complex :math:`2\times 2` matrix :math:`Y_{km}`, the admittance-matrix.  We also have a value :math:`L_{km}` (the branch limit).

2. For every bus :math:`k` we have two positive values, :math:`L_k \, \text{and} \, U_k`, the voltage limits.

3. For every bus :math:`k` we have two values, :math:`P^{d}_k \, \text{and} \, Q^{d}_k`, the active and reactive loads (demands).

4. For every generator :math:`i` we have the active and reactive generation limits, :math:`P^{\min}_i, \, P^{\max}_i, \, Q^{\min}_i, \, \text{and} \, Q^{\max}_i`, as well as a (convex) function :math:`F(i)`.

5. At a bus :math:`k` we have a set :math:`G(k)`, the set of generators located at the bus.
6. For a bus :math:`k`, denote by :math:`\delta^+(k) \, \text{and} \, \delta^-(k)` the set of branches of the form :math:`km \, \text{and} \, mk`, respectively.

7. Other noteworthy details: there may be multiple parallel branches :math:`km` but in this writeup, for the sake of simplicity, this is being ignored (the mathematical description is easily adapted).  There may be multiple generators at a given bus (we account for this explicitly).  The same bus may generators and nonzero load.  In this simple description we ignore bus shunts -- the mathematical description is also easily modified to account for these.

The following variables are used in this problem.

1. For each bus :math:`k` we have a complex voltage :math:`V_k = |V_k| e^{j \theta_k}` (where :math:`j = \sqrt{-1}`).

2. For each branch :math:`km` we have the two complex variables :math:`S_{km} \, \text{and} \, S_{mk}`, the complex power injected into the branch at :math:`k \, \text{and at} \, m`, respectively. (Recall that branches are given as ordered
   pairs).

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


The objective to optimize is the following:



.. math::
   \min \sum_{\text{generators} \, i} F_i(P^g_i).

.. testsetup:: mod
    # Set pandas options for displaying dataframes, if needed
    import pandas as pd
    pd.options.display.max_rows = 10

.. tabs::

    .. tab:: ``availability``

        Give interpretation of input data.

        .. doctest:: mod
            :options: +NORMALIZE_WHITESPACE

            >>> from gurobi_optimods import datasets
            >>> data = datasets.load_workforce()
            >>> data.availability
               Worker      Shift
            0     Amy 2022-07-02
            1     Amy 2022-07-03
            2     Amy 2022-07-05
            3     Amy 2022-07-07
            4     Amy 2022-07-09
            ..    ...        ...
            67     Gu 2022-07-10
            68     Gu 2022-07-11
            69     Gu 2022-07-12
            70     Gu 2022-07-13
            71     Gu 2022-07-14
            <BLANKLINE>
            [72 rows x 2 columns]

        In the model, this corresponds to ...

    .. tab:: ``shift_requirements``

        Another bit of input data (perhaps a secondary table)

|

Code
----

Self contained code example to run the mod from an example dataset. Example
datasets should bd included in the ``gurobi_optimods.datasets`` module for
easy access by users.

.. testcode:: mod

    from gurobi_optimods.opf import solve_opf_model, read_settings_from_file, read_case_from_file
    from gurobi_optimods.datasets import load_caseopf

    settings = {"doac": True, "use_ef": True}
    # load path to case file
    casefile = load_caseopf("9")
    # read case file and return a case dictionary
    case = read_case_from_file(casefile)
    # solve opf model and return a solution and the final objective value
    solution, objval = solve_opf_model(settings, case)

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
