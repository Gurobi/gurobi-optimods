DCOPF Problem Formulation
=========================

Problem Specification
---------------------


.. tabs::

   .. tab:: Domain-Specific Description

        In the so-called DC approximation to the standard ACOPF problem it is assumed that all voltage magnitudes are equal to :math:`1.0` and that across all branches the phase angle difference is very small.  The active power flow equations are linearized, using these assumptions, and the reactive power flow constraints are ignored.  The objective function is the same as for ACOPF. In summary we obtain a linear approximation (not a relaxation) to standard ACOPF which is very commonly used in energy markets.

	The network data requirements are likewise limited.  For each branch :math:`km` the models require the branch reactance :math:`x_{km}` as well as a ratio :math:`\tau_{km}` and angle :math:`\phi_{km}`; the latter two are relevant only in the case of transformers.  In the non-transformer case we assume :math:`\tau_{km} = 1` and :math:`\phi_{km} = 0`.

   .. tab:: Optimization Model
       The following variables are used in the DCOPF problem.

       #. For each bus :math:`k` we have a phase angle :math:`\theta_k`.

       #. For each branch :math:`km` we have the real power flow :math:`P_{km}` from :math:`k` to :math:`m`. As we shall see from the model, there is no need for a corresponding variable :math:`P_{mk}`.

       #. For each generator :math:`i` there is the (active) generation :math:`P^{g}_i`.

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

   .. tab:: Branch-switching

       An important binary MILP extension of DCOPF is that where a binary variable :math:`z_{km}` is used to decide if a branch :math:`km` is "on" (:math:`z_{km} = 1`) or not (:math:`z_{km} = 0`).  To achive this goal, we simply reformulate the power flow definition as

       .. math::

        P_{km} = \frac{z_{km}}{\tau_{km} \,x_{km}}(\theta_k - \theta_m - \phi_{km}).

       This constrained is linearized through standard methods.
