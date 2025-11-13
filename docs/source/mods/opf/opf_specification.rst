Full Problem Specification
==========================


As input data for an OPF problem we have a power system (*power grid* for non-engineers) consisting of a network of buses (nodes) and branches (arcs).

Each branch is given as an ordered pair :math:`km`, where :math:`k` is the *from* bus and :math:`m` is the *to* bus.

1. For each branch :math:`km` we have a complex :math:`2\times 2` matrix :math:`Y_{km}`, the admittance-matrix. The following notation will be used below

       .. math::

        Y _{km} =  \left( \begin{array}{c c}
        G_{kk} +j B_{kk} & G_{km} +j B_{km}\\
        G_{mk} +j B_{mk} & G_{mm} +j B_{mm}
        \end{array} \right).

   We also have a value :math:`L_{km}` defining the branch limit.

2. For every bus :math:`k` we have two positive values, :math:`L_k \, \text{and} \, U_k` with :math:`L_k \leq U_k`, defining the voltage limits.

3. For every bus :math:`k` we have two values, :math:`P^{d}_k \, \text{and} \, Q^{d}_k`, the active and reactive loads (demands), respectively.

4. For every generator :math:`i` we have the active and reactive generation limits, :math:`P^{\min}_i, \, P^{\max}_i, \, Q^{\min}_i, \, \text{and} \, Q^{\max}_i` with :math:`P^{\min}_i, \leq P^{\max}_i, \, Q^{\min}_i, \leq Q^{\max}_i`, as well as a (convex) function :math:`F_i` defining the generation cost.

5. At a bus :math:`k` we have a set :math:`G(k)` defining the set of generators located at bus :math:`k`.

6. For a bus :math:`k`, denote by :math:`\delta^+(k) \, \text{and} \, \delta^-(k)` the set of branches outgoing and incoming branches of the form :math:`km \, \text{and} \, mk`, respectively.

7. There may be multiple parallel branches :math:`km` but in this writeup, for the sake of simplicity, this is being ignored (the mathematical description is easily adapted).  There may be multiple generators at a given bus (we account for this explicitly).  The same bus may have generators and a nonzero load (demand).  In this simple description we ignore bus shunts.

8. Associated with each branch :math:`km` there are values :math:`\tau_{km} > 0` (ratio) and :math:`\phi_{km}` (angle).  These are meaningful in the case of transformers; for a non-transformer branch we set :math:`\tau_{km} = 1` and :math:`\phi_{km} = 0`.

9. In modern versions discussed in the :ref:`Recommended Literature <recommended-label>`, network devices such as FACTS, phase-shifters, and impedance corrections are also incorporated into the model. We do not account for these in this mod.


.. _acopf-label:

ACOPF
-----

In the standard ACOPF problem we set a minimum-cost operating point (complex voltages at all buses, active and reactive power generation at each generator) for a power system so that the correct amount of active and reactive power is delivered to each bus and all branch limits are satisfied; both using the AC power flow laws.  Under this model, cost is incurred at the generators.  At each generator we are given a convex quadratic or piecewise-linear (convex) cost function associated with active power generation.


AC Optimization Model
~~~~~~~~~~~~~~~~~~~~~

The following variables are used in the ACOPF problem.

    - For each bus :math:`k` we have a complex voltage :math:`V_k = |V_k| e^{j \theta_k}` (where :math:`j = \sqrt{-1}`).

    - For each branch :math:`km` we have the two complex variables :math:`S_{km} \, \text{and} \, S_{mk}`, the complex power injected into the branch at :math:`k \, \text{and at} \, m`, respectively. Recall that branches are given as ordered pairs.

    - For each generator :math:`i` there is a (complex) generation :math:`P^{g}_i + j Q^{g}_i`.

There are five classes of constraints:

    #. branch AC power definitions
    #. flow balance
    #. branch limits
    #. generator limits
    #. voltage bounds

1. The first set of constraints embody AC power flow laws. For each branch :math:`km` we have

.. math::

    S _{km} =  V_k \left[ Y \left( \begin{array}{r}
    V_{k} \\
    V_{m}
    \end{array} \right) \right]^*_k \quad \text{and} \quad S _{mk} =  V_m \left[ Y \left( \begin{array}{r}
    V_{k} \\
    V_{m}
    \end{array} \right) \right]^*_m,

where * indicates the complex conjugate.  These are the only nonconvex expressions in the entire model.  We will also describe alternative representations of these constraints using different notation.

2. The next set of constraints represent flow balance.  For each branch :math:`k` we have

.. math::

    \sum_{i \in G(k)} (P^{g}_i + jQ^{g}_i) \ - \ (P^{d}_k + j Q^{d}_k) \quad = \quad \sum_{km \in \delta^+(k)} S_{km} \ + \ \sum_{mk \in \delta^-(k)}S_{km}.


In this expression, the left-hand side is the net excess of generation over load at bus :math:`k`.  The right-hand side is the total power injected into the grid at bus :math:`k`.

3. Next we have branch limits.  For each branch :math:`km`,

.. math::

    \max\{ |S_{km}|, |S_{mk}| \} \ \le \, L_{km}.

Note: there are alternative versions of this particular constraint.  In a power system under normal (non-stressed) operations this constraint is slack for all but a small number of branches.

4.-5. Finally we have generator limits and voltage bounds.

.. math::

    P^{\min}_i \, \le \, P^{g}_i \, \le \, P^{\max}_i \quad \text{and} \quad  Q^{\min}_i \, \le \, Q^{g}_i \, \le \, Q^{\max}_i, \ \text{for each generator} \, i,

and

.. math::
    L_k \, \le \, |V_k| \, \le \, U_k, \ \text{for each bus} \, k.


The objective to optimize is the following cost expression:


.. math::
    \min \sum_{\text{generators} \, i} F_i(P^g_i).


Above we provided an expression defining the complex power flow :math:`S_{km}` in terms of the admittance matrix :math:`Y`. In the following we derive a formulation independent of complex numbers. Denote

.. math::
    \theta_{km} = \theta_k - \theta_m - \phi_{km},

The real and imaginary components of :math:`S_{km}` are

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

the power flow quantities can be rewritten without the need of complex numbers as

.. math::
    P_{km} \ = \ G_{kk} v_k^{(2)} + G_{km} c_{km} + B_{km} s_{km}

and

.. math::
    Q_{km} \ = \ -B_{kk} v_k^{(2)}  - B_{km} c_{km} + G_{km} s_{km}.

respectively.

.. _polar-label:

Polar Formulation
~~~~~~~~~~~~~~~~~

The :ref:`ACOPF formulation<acopf-label>` described above using :math:`\cos(\theta_{km})` and :math:`\sin(\theta_{km})` terms is called the *polar* formulation. Despite its high nonlinearity, it is the most widely used formulation for solving ACOPF models. This is the default formulation that Gurobi solves in this mod.

.. _qcqp-label:

QCQP
~~~~

ACOPF can be reformulated as a Quadratically Constraints Quadratic Program (QCQP) by performing two reformulation steps.

1. For each bus :math:`k`, introduce the real variables :math:`e_k \, \text{and} \, f_k` and set :math:`v^{(2)}_k \, = \, e_k^2 + f_k^2`.

2. For each branch :math:`km`, set :math:`c_{km} = e_k e_m + f_k f_m` and :math:`s_{km} = -e_k f_m + f_k e_m`.

These constraints render an exact reformulation rendering the problem as a QCQP, i.e., one that removes the sines and cosines from the formulation.  We remind the reader that, for example, we are writing the active power flow injected at bus :math:`k` on branch :math:`km` through the constraints

.. math::
    P_{km} \ = \ G_{kk} v_k^{(2)} + G_{km} c_{km} + B_{km} s_{km},

(and similarly with :math:`Q_{km}, \, P_{mk}, \,\text{and} \,Q_{mk}`).

This is the so-called *cartesian* (or *rectangular*) formulation for ACOPF.

.. _jabr-label:

Jabr Relaxation
~~~~~~~~~~~~~~~

We can obtain a Second-Order Cone (SOC) relaxation of the ACOPF formulation
by introducing auxiliary variables :math:`v^{(2)}_k, \ c_{km} \ \text{and} \ s_{km}`, removing the nonconvex definitions of such variables (which involve cosines and sines) and adding the rotated cone constraints

.. math::
    c_{km}^2 \ + \ s_{km}^2 \ \le \ v_k^{(2)} v_m^{(2)} \ \text{for every branch} \, km.


The resulting relaxation can prove very tight, though, despite its convexity, challenging in large cases.

.. _dcopf-label:

DCOPF
-----

In the so-called DC approximation to the standard ACOPF problem it is assumed that all voltage magnitudes are equal to :math:`1.0` and that across all branches the phase angle difference is very small.  The active power flow equations are linearized, using these assumptions, and the reactive power flow constraints are ignored.  The objective function is the same as for ACOPF. In summary we obtain a linear approximation (not a relaxation) to standard ACOPF which is very commonly used in energy markets.

The network data requirements are likewise limited.  For each branch :math:`km` the models require the branch reactance :math:`x_{km}` as well as a ratio :math:`\tau_{km}` and angle :math:`\phi_{km}`; the latter two are relevant only in the case of transformers.  In the non-transformer case we assume :math:`\tau_{km} = 1` and :math:`\phi_{km} = 0`.


DC Optimization Model
~~~~~~~~~~~~~~~~~~~~~

The following variables are used in the DCOPF problem.

    - For each bus :math:`k` we have a phase angle :math:`\theta_k`.

    - For each branch :math:`km` we have the real power flow :math:`P_{km}` from :math:`k` to :math:`m`. As we shall see from the model, there is no need for a corresponding variable :math:`P_{mk}`.

    - For each generator :math:`i` there is the (active) generation :math:`P^{g}_i`.

There are four classes of constraints:

    #. branch power definitions
    #. flow balance
    #. branch limits
    #. generator limits

1. The first set of constraints embodies DC power flow laws. For each branch :math:`km` we have

.. math::

    P_{km} = \frac{1}{\tau_{km} \,x_{km}}(\theta_k - \theta_m - \phi_{km}).

It is (implicitly) assumed that :math:`P_{mk} = -P_{km}`.

2. The next set of constraints represents active flow balance at each bus :math:`k`.

.. math::

    \sum_{i \in G(k)} P^{g}_i \ - \ P^{d}_k \quad = \quad \sum_{km \in \delta^+(k)} P_{km} \ - \ \sum_{mk \in \delta^-(k)}P_{mk}.


In this expression, the left-hand side is the net excess of generation over load at bus :math:`k`.  The right-hand side is the total power injected into the grid at bus :math:`k`.

3. Next we have branch limits.  For each branch :math:`km`,

.. math::

    |P_{km}| \ \le \, L_{km}.

Note: there are alternative versions of this particular constraint.  In a power system under normal (non-stressed) operations this constraint is slack for all but a small number of branches.

4. Finally we have generator limits.

.. math::

    P^{\min}_i \, \le \, P^{g}_i \, \le \, P^{\max}_i, \ \text{for each generator} \, i.

The objective to optimize is the same (usually quadratic convex) cost expression as for ACOPF:


.. math::
    \min \sum_{\text{generators} \, i} F_i(P^g_i).

.. _branchswitching-label:

Branch Switching
~~~~~~~~~~~~~~~~

An important binary MILP extension of DCOPF is that where a binary variable :math:`z_{km}` is used to decide if a branch :math:`km` is "on" (:math:`z_{km} = 1`) or not (:math:`z_{km} = 0`).  To achive this goal, we simply reformulate the power flow definition as

.. math::

    P_{km} = \frac{z_{km}}{\tau_{km} \,x_{km}}(\theta_k - \theta_m - \phi_{km}).

This constrained is linearized through standard methods.

Branch switching is most often used in DCOPF only due to its additional complexity. Still, it can be applied to an ACOPF model as well. Please note that introducing branch-switching greatly increases problem complexity and thus, obtaining feasible solution to branch-switching ACOPF models is very difficult.
