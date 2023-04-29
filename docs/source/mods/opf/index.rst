Optimal Power Flow
==================

The operation of power systems relies on a number of optimization tasks, known as 'optimal power flow' of OPF problems.  Here we assume basic familiarity with concepts such as 'voltage' (potential energy), 'current' (charge flow) and 'power' (instantaneous energy generation or consumption).  The engineering community also
uses the term 'bus' (nodes in a network, with some simplification) and 'branch' (a connection between two buses, typically a line or a transformer).

To learn more
---------------------

- :footcite:t:`andersson1`
      - Thorough survey of relaxations for power flow problems.

- :footcite:t:`bergenvittal`
      - Thorough survey of relaxations for power flow problems.

- :footcite:t:`bienstockpower`
      - Introduction to various power related problems with emphasis on modeling and optimization.

- :footcite:t:`molzahnhiskens`
      - Thorough survey of relaxations for power flow problems.

- :footcite:t:`overbyebook`
      - Textbook on power systems.


Problem Specification
---------------------

As input data for an OPF problem we have a power system ('power grid' for non-engineers) consisting of a network of buses and branches. Each branch is given as
an ordered pair :math:`km`, where :math:`k` is the "from" bus and :math:`m` is the "to" bus.

1. For each branch :math:`km` we have a complex :math:`2\times 2` matrix :math:`Y_{km}`, the admittance-matrix. The following notation will be used below

       .. math::

           Y _{km} =  \left( \begin{array}{c c}
	   G_{kk} +j B_{kk} & G_{km} +j B_{km}\\
	   G_{mk} +j B_{mk} & G_{mm} +j B_{mm}
	   \end{array} \right).

   We also have a value :math:`L_{km}` (the branch limit).

2. For every bus :math:`k` we have two positive values, :math:`L_k \, \text{and} \, U_k`, the voltage limits.

3. For every bus :math:`k` we have two values, :math:`P^{d}_k \, \text{and} \, Q^{d}_k`, the active and reactive loads (demands).

4. For every generator :math:`i` we have the active and reactive generation limits, :math:`P^{\min}_i, \, P^{\max}_i, \, Q^{\min}_i, \, \text{and} \, Q^{\max}_i`, as well as a (convex) function :math:`F(i)`.

5. At a bus :math:`k` we have a set :math:`G(k)`, the set of generators located at the bus.
6. For a bus :math:`k`, denote by :math:`\delta^+(k) \, \text{and} \, \delta^-(k)` the set of branches of the form :math:`km \, \text{and} \, mk`, respectively.

7. Other noteworthy details: there may be multiple parallel branches :math:`km` but in this writeup, for the sake of simplicity, this is being ignored (the mathematical description is easily adapted).  There may be multiple generators at a given bus (we account for this explicitly).  The same bus may generators and nonzero load.  In this simple description we ignore bus shunts -- the mathematical description is also easily modified to account for these.

8. Associated with each branch :math:`km` there are values :math:`\tau_{km} > 0` and :math:`\phi_{km}`.  These are meaningful in the case of transformers; for a non-transformer branch we write :math:`\tau_{km} = 1` and :math:`\phi_{km} = 0`.

9. In modern versions described below, network devices (FACTS, phase-shifters, impedance corrections are also incorporated into the model.



More complete descriptions
--------------------------

.. toctree::
   :maxdepth: 1

   acopf
   dcopf

Basic use examples
------------------

.. testsetup:: mod

    # Set pandas options for displaying dataframes, if needed
    import pandas as pd
    pd.options.display.max_rows = 10

.. tabs::

    .. tab:: ``Reading a case``

        Here we show how to read in a 'case' based on the NY State network.

        .. doctest:: mod
            :options: +NORMALIZE_WHITESPACE

	       >>> import gurobipy as gp
	       >>> from gurobi_optimods.opf import read_case_from_mat_file
	       # This is a basic module used to read a casefile
	       >>> casefile = 'caseNY.mat' #make sure you have the right path
	       >>> case = read_case_from_mat_file(casefile)  #creates main dictionary

	       #example of a bus:
	       >>> buses = case['bus']
	       >>> buses[101]
	       {'bus_i': 3600, 'type': 2, 'Pd': 2.14, 'Qd': 1.16, 'Gs': 0.0, 'Bs': 0.0, 'area': 24.0, 'Vm': 1.0163928, 'Va': -91.575965, 'baseKV': 138.0, 'zone': 1.0, 'Vmax': 1.1, 'Vmin': 0.9}


	       #example of a branch:
	       >>> branches = case['branch']
	       >>> branches[10]
	       {'fbus': 3362.0, 'tbus': 3542.0, 'r': 0.000676, 'x': 0.012043, 'b': 0.12494, 'rateA': 1046.0, 'rateB': 1098.3, 'rateC': 1150.6, 'ratio': 1.0, 'angle': 0.0, 'status': 1.0, 'angmin': 0.0, 'angmax': 0.0}

	       #example of a generator
	       >>> generators = case['gen']
	       >>> gens[5]
	       {'bus': 3505, 'Pg': 487.92, 'Qg': 0.0, 'Qmax': 0.0, 'Qmin': 0.0, 'Vg': 1.02, 'mBase': 100.0, 'status': 1, 'Pmax': 883.0, 'Pmin': 1.0, 'Pc1': 0.0, 'Pc2': 0.0, 'Qc1min': 0.0, 'Qc1max': 0.0, 'Qc2min': 0.0, 'Qc2max': 0.0, 'ramp_agc': 0.0, 'ramp_10': 0.0, 'ramp_30': 0.0, 'ramp_q': 0.0, 'apf': 883.0}

.. footbibliography::
