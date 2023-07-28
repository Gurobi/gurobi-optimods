Optimal Power Flow
==================

.. toctree::
    :hidden:

    opf_specification

The operation of power systems relies on a number of optimization tasks, known
as Optimal Power Flow (OPF) problems. The objective of a standard OPF problem is
to minimize operational cost such that the underlying grid constraints on
generation, demand, and voltage limits, are satisfied.

In this mod, we consider the cases of Alternating Current (AC) and Direct
Current (DC) OPF formulations. The ACOPF problem in its natural form requires
the introduction of complex numbers in order to formulate voltage constraints.
This Mod uses a cartesian-coordinates formulation of ACOPF that reformulates
complex-valued terms via nonconvex quadratic relationships. The DCOPF problem is
an approximation of the ACOPF problem where additional assumptions to produce
linear constraints. While the additional assumptions result in potential loss of
solution accuracy, they make the DCOPF problem a lot easier to solve. This is
especially useful if the solution accuracy can be neglected in favor of solution
time and problem size. For full details of the formulations, please refer to the
:doc:`opf_specification`.

Here we assume basic familiarity with concepts such as *voltage* (potential
energy), *current* (charge flow), and *power* (instantaneous energy generation
or consumption). The engineering community also uses the terms *bus* to refer to
nodes in a network, and branch to refer to arcs in a network (a connection
between two buses, typically a line or a transformer). For more details and
comprehensive descriptions of power systems and the underlying problems, please
refer to the `Recommended Literature`_ section.

Solving an OPF Problem
----------------------

This Mod has multiple API functions. Each function takes a dictionary as input
which describes an OPF case. This case dictionary holds all essential
information about the underlying network: buses, branch connections, and
generators. The case dictionary follows the `MATPOWER Case Format conventions
<https://matpower.org/docs/ref/matpower7.1/lib/caseformat.html>`_. The details
of the input format are discussed in the `Case and Result Dictionaries`_
section.

Several pre-defined MATPOWER cases can be read in from the
`gurobi_optimods.datasets` module. The following code loads a 9 bus grid
example:

.. testcode:: opf

    from gurobi_optimods import datasets

    case = datasets.load_opf_example("case9")

After reading in or otherwise generating a case dictionary, we can solve an OPF
problem defined by the given network data. For this task, we use the
:func:`gurobi_optimods.opf.solve_opf` function. We can define the type of the
OPF problem that we want to solve by defining the ``opftype`` argument when
calling the function. Currently, the available options are ``AC``, ``AC_relax``,
and ``DC``.

- The ``AC`` setting solves an ACOPF problem defined by the given network data.
  The ACOPF problem is formulated as a nonconvex bilinear model as described in
  the :ref:`ACOPF <acopf-label>` section of the :doc:`opf_specification`. This
  setting yields the most accurate model of the physical power system. However,
  it is also the most difficult problem to solve and thus usually leads to the
  longer runtimes.

- The ``AC_relax`` setting solves a Second Order Cone (SOC) relaxation of the
  nonconvex bilinear ACOPF problem formulation defined by the given network
  data. The relaxation is constructed by dropping nonconvex bilinear terms but
  simultaneously keeping the convex JABR inequalities, see :ref:`JABR Relaxation
  <jabr-label>` for more details. This setting often yields a good approximation
  of the physical power system and is of moderate difficulty.

- The ``DC`` setting solves a DCOPF problem defined by the given network data.
  The DCOPF problem is a linear approximation of the ACOPF problem, see
  :ref:`DCOPF <dcopf-label>` section of the :doc:`opf_specification` for more
  details. This setting only yields a crude approximation of the physical power
  system, but is usually an easy problem that can be solved very quickly even
  for large networks.

The ``solve_opf`` function solves an ``AC`` problem unless the ``opftype``
argument specifies otherwise.

.. testcode:: opf

    from gurobi_optimods import opf

    result = opf.solve_opf(case, opftype="AC")

.. testoutput:: opf
    :options: +NORMALIZE_WHITESPACE +ELLIPSIS

    ...
    Optimize a model with 73 rows, 107 columns and 208 nonzeros
    ...
    Optimal solution found...
    ...
    Objective value = 5296...
    ...

ACOPF and `Branch-Switching`_ models are most often very hard to solve to
optimality. For this reason, it is best to pass specific solver settings such
as, e.g., a `TimeLimit
<https://www.gurobi.com/documentation/current/refman/timelimit.html>`_. This can
be done by using the ``solver_params`` argument. For a full list of all Gurobi
parameters please refer to `our documentation
<https://www.gurobi.com/documentation/current/refman/parameter_descriptions.html>`_.

.. testcode:: opf

    result = opf.solve_opf(
        case,
        opftype="AC",
        solver_params={"TimeLimit": 60}
    )

.. testoutput:: opf
    :hide:
    :options: +NORMALIZE_WHITESPACE +ELLIPSIS

    Set parameter TimeLimit to value 60
    ...
    Optimal solution found...
    ...

The Mod returns the result as a dictionary, following the same `MATPOWER Case
Format conventions
<https://matpower.org/docs/ref/matpower7.1/lib/caseformat.html>`_ as the case
dictionary. However, in the result dictionary some object entries are modified
compared to the input case dictionary. These modified fields hold the solution
values of the optimization. In some cases, there are also additional fields to
store the solution information. We discuss all details of the result dictionary
in the `Case and Result Dictionaries`_ section.

Graphical Representation of Feasible Solutions
----------------------------------------------

In addition to solving an OPF problem, this Mod also provides plotting functions
to display graphical representation of the network and the OPF result. There are
already very involved graphical tools to represent OPF solutions provided by
other packages such as:

- `MATPOWER <https://matpower.org>`_
- `PyPSA <https://pypsa.org/>`_
- `pandapower <http://www.pandapower.org/>`_

thus the graphical representation provided by this Mod is kept intentionally
simple. In order to use this functionality, it is necessary to install the
``plotly`` package as follows::

    pip install plotly


Plotting the Result
~~~~~~~~~~~~~~~~~~~

In order to plot a previously obtained result, you must provide :math:`(x, y)`
coordinates for all buses in the network. Coordinates are provided as a
dictionary mapping bus IDs to coordinates. The OptiMods datasets module provides
an example set of coordinates for plotting the 9 bus test case:

.. TODO is it mapping bus IDs? Or positions in the bus list of the case?

.. doctest:: opf
    :options: +NORMALIZE_WHITESPACE

    >>> coordinates = datasets.load_opf_extra("case9-coordinates")
    >>> from pprint import pprint
    >>> pprint(coordinates)
    {1: (44.492, -73.208),
     2: (41.271, -73.953),
     3: (41.574, -73.966),
     4: (40.814, -72.94),
     5: (43.495, -76.451),
     6: (42.779, -78.427),
     7: (44.713, -73.456),
     8: (43.066, -76.214),
     9: (43.048, -78.854)}

Given a solution and a coordinate mapping, the plotting functions return plotly
figures which can be displayed in a web browser. In the following code, we solve
the DCOPF problem for a real-world dataset for the city of New York and produce
a solution plot::

    case = datasets.load_opf_example("caseNY")
    solution = opf.solve_opf(case, opftype='DC')
    coordinates = datasets.load_opf_extra("caseNY-coordinates")
    fig = opf.solution_plot(case, coordinates, solution)
    fig.show()  # open plot in a browser window

.. NB: we cannot do this in a doctest as the model is large and requires a
.. license. So it's tricky to ensure this code stays working ...

.. figure:: ../figures/opf.png

    DCOPF solution for the New York power grid example dataset

The above image shows the grid solution generated from the given network data,
plotted using the provided coordinates. The colored circles depict generators
and the amount of power they generate:

- Black bus: Power generation :math:`\leq 75` and load :math:`< 50`
- Blue bus: Power generation :math:`\leq 75` and load :math:`\geq 50`
- Purple bus: Power generation :math:`> 75`
- Orange bus: Power generation :math:`> 150`
- Red bus: Power generation :math:`> 500`

Branch-Switching
----------------

An important extension of the OPF problem is Branch Switching, where branches
may be turned off. Note that already turning off a single branch changes the
whole power flow through the network. Thus in practice, it is rare that branches
are turned off at all, even if this option is enabled. If any are turned off,
then it is usually only a small fraction of the overall power grid. For the
mathematical formulation, please refer to the :ref:`Branch-Switching
<branchswitching-label>` subsection of the :doc:`opf_specification`.

To enable branch-switching in a given OPF problem, set the ``branch_switching``
argument to ``True`` when calling :func:`gurobi_optimods.opf.solve_opf`. The Mod
additionally offers the possibility to control the number of branches that must
remain switched on via the ``min_active_branches`` argument. In practice, it is
expected that only a very small fraction of branches are turned off. Thus, the
default value of the ``min_active_branches`` argument is 0.9 (90%). In the
following example, we solve a modified version of the 9 bus network to see
whether branch switching allows a better solution.

.. testcode:: opf

    case = datasets.load_opf_example("case9-switching")
    result = opf.solve_opf(
        case, opftype="AC", branch_switching=True, min_active_branches=0.1
    )

.. testoutput:: opf
    :hide:
    :options: +NORMALIZE_WHITESPACE +ELLIPSIS

    ...
    Optimize a model with 278 rows, 185 columns and 694 nonzeros
    ...

Plotting the resulting solution shows that one branch has been turned off in the
optimal solution. Please note, that the used examplary network has been
artificially adjusted to achieve this result and this is **not** the usual
behavior in a realistic power grid of such small size.

.. figure:: ../figures/switching_opf.png

    Branch switching solution. The plot highlighted the switched-off branch in
    red

Violations for Pre-defined Voltage Values
-----------------------------------------

In practice we may have voltage magnitudes and voltage angles for each bus at
hand and would like to know whether these values are actually feasible within a
given network. To tackle this problem, we can use the
:func:`gurobi_optimods.opf.compute_violations` function. This function takes a
set of bus voltages in addition to the case data and returns the computed
violations in this voltage solution.

Bus voltages are provided as a dictionary mapping the Bus ID to a pair
:math:`(V_m, V_a)` where :math:`V_m` is the voltage magnitude and :math:`V_a` is
the voltage angle. An example is provided for the 9 bus case:

.. doctest:: opf

    >>> voltages = datasets.load_opf_extra("case9-voltages")
    >>> pprint(voltages)
    {1: (1.089026, 0.0),
     2: (1.099999, 20.552543),
     3: (1.090717, 16.594399),
     4: (1.084884, -2.408447),
     5: (1.096711, 2.43001),
     6: (1.099999, 11.859651),
     7: (1.072964, 9.257936),
     8: (1.066651, 11.200108),
     9: (1.08914, 2.847507)}


Using this voltage data, we can check for possible model violations by calling
the :func:`gurobi_optimods.opf.compute_violations` function. The function
returns a dictionary which follows the `MATPOWER Case Format
<https://matpower.org/docs/ref/matpower7.1/lib/caseformat.html>`_ with
additional fields storing the violations for particular buses and branches.

The following fields in the violations dictionary are added to store violations
data:

- ``violation['bus'][i]['Vmviol']`` Voltage magnitude violation at bus `i`
- ``violation['bus'][i]['Pviol']`` real power injection violation at bus `i`
- ``violation['bus'][i]['Qviol']`` reactive power injection violation at bus `i`
- ``violation['branch'][i]['limitviol']`` branch limit violation at branch `i`

.. testcode:: opf

    volts_dict = datasets.load_opf_extra("case9-voltages")
    case = datasets.load_opf_example("case9")
    violations = opf.compute_violations(case, volts_dict)

.. testoutput:: opf
    :hide:
    :options: +NORMALIZE_WHITESPACE +ELLIPSIS

    ...
    Checking flow balance constraints.
    ...

.. doctest:: opf

    >>> print(violations['branch'][6]['limitviol'])
    66.33435016796234
    >>> print(violations['bus'][3]['Pviol'])
    -318.8997836192236

In this case, the limit at branch 6 and the real power injection at bus 3 are
violated by the given input voltages.

Inspecting Violations Graphically
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Similar to generating a graphical representation of a feasible solution, it is also possible to generate a figure representing the violations within a given power grid. We can use the :func:`gurobi_optimods.opf.violation_plot` to generate a :class:`plotly.graph_objects.Figure` object that can be displayed in, e.g., a browser window. Please note that it is required to install the ``plotly`` package to use this functionality. In addition to bus coordinates and case information, we also need to provide the violations dictionary to the function. In the following we use the previously discussed violated solution and plot the result.

.. code-block::

    volts_dict = datasets.load_opf_extra("case9-voltages")
    case = datasets.load_opf_example("case9")
    coords_dict = datasets.load_opf_extra("case9-coordinates")
    violations = opf.compute_violations(case, volts_dict)
    fig = opf.violation_plot(case, coords, violations)
    fig.show()

.. image:: ../figures/violations_opf.png

In the above image, you can see the power grid generated out of the given network data together with the coordinate and violation information. The red circles depict buses where the voltage magnitude or real or reactive power injections are violated. Red marked branches depict branches with violated limits.


Case and Result Dictionaries
----------------------------

This mod uses so-called *case* and *result* dictionaries for input and output. Both dictionaries are meant to follow the `MATPOWER Case Format <https://matpower.org/docs/ref/matpower7.1/lib/caseformat.html>`_ and are described in full detail below.

Case Dictionary
~~~~~~~~~~~~~~~

The case dictionary for this mod expects a dictionary with keys ``baseMVA``, ``bus``, ``branch``, ``gen``, and ``gencost``. All other entries of the dictionary are ignored in the current version of the mod. The value stored via the ``baseMVA`` key is a numerical float value. The values stored in the case dictionary via keys ``bus``, ``branch``, ``gen``, and ``gencost`` are lists of dictionaries, where each dictionary holds specific data about the particular object. Every single object is defined by a dictionary holding entries following the `MATPOWER Case Format <https://matpower.org/docs/ref/matpower7.1/lib/caseformat.html>`_, e.g., every bus has a bus number ``bus_i``, real power demand ``Pd``, etc.

.. testcode:: opf

    case = datasets.load_opf_example("case9")

.. testoutput:: opf
    :hide:
    :options: +NORMALIZE_WHITESPACE

    ...

.. doctest:: opf
    :options: +NORMALIZE_WHITESPACE

    >>> case['baseMVA']
    100.0

    >>> buses = case['bus']
    >>> buses[0]
    {'bus_i': 1, 'type': 3, 'Pd': 0.0, 'Qd': 0.0, 'Gs': 0.0, 'Bs': 0.0, 'area': 1.0, 'Vm': 1.0, 'Va': 0.0, 'baseKV': 345.0, 'zone': 1.0, 'Vmax': 1.1, 'Vmin': 0.9}

    >>> branches = case['branch']
    >>> branches[1]
    {'fbus': 4, 'tbus': 5, 'r': 0.017, 'x': 0.092, 'b': 0.158, 'rateA': 250.0, 'rateB': 250.0, 'rateC': 250.0, 'ratio': 0.0, 'angle': 0.0, 'status': 1.0, 'angmin': -360.0, 'angmax': 360.0}

    >>> generators = case['gen']
    >>> generators[2]
    {'bus': 3, 'Pg': 85, 'Qg': 0, 'Qmax': 300, 'Qmin': -300, 'Vg': 1, 'mBase': 100, 'status': 1, 'Pmax': 270, 'Pmin': 10, 'Pc1': 0, 'Pc2': 0, 'Qc1min': 0, 'Qc1max': 0, 'Qc2min': 0, 'Qc2max': 0, 'ramp_agc': 0, 'ramp_10': 0, 'ramp_30': 0, 'ramp_q': 0, 'apf': 0}

    >>> generatorcosts = case['gencost']
    >>> generatorcosts[2]
    {'costtype': 2.0, 'startup': 3000.0, 'shutdown': 0.0, 'n': 3.0, 'costvector': [0.1225, 1.0, 335.0]}

There is also the convenience function :func:`gurobi_optimods.opf.read_case_matpower` which reads in a standard MATLAB ``.mat`` data file holding the network data. The data stored in the ``.mat`` file has to follow the `MATPOWER Case Format conventions <https://matpower.org/docs/ref/matpower7.1/lib/caseformat.html>`_ in order to be accepted by the function. This function returns a case dictionary which can be read by the solver methods::

    case = opf.read_case_matpower("my_case.mat")
    solution = solve_opf(case, opftype="AC")


.. warning::

    In the current version of the mod, we only accept generator costs with `costtype = 2`, i.e., polynomial model, up to degree 2, i.e, `n=3` in the `gencost` structure. For now, these seem to be the most commonly used settings in practice. If a different costtype or `n` value is provided, an error is issued. It is possible that more functionality will be added in a future release.


Result Dictionary
~~~~~~~~~~~~~~~~~

The following fields in the result dictionary are altered or added to store solution data.

- The field ``result['et']`` holds the runtime value of the whole solution process in seconds.
- The field ``result['success']`` defines whether at least one feasible solution has been found (1) or no feasible solution is available (0).
- The field ``result['f']`` holds the solution objective value (only valid if ``result['success'] == 1``).
- If a feasible solution has been found, the ``bus`` entries

  * ``result['bus'][i]['Vm']``
  * ``result['bus'][i]['Va']``

  store the voltage magnitude (Vm) and voltange angle (Va) values in the optimal solution for bus `i`.
- If we solved a DCOPF problem, additional fields ``result['bus'][i]['mu']`` hold the shadow prices for balance constraints at bus `i`.
- If a feasible solution has been found, the ``gen`` entries

  * ``result['gen'][i]['Pg']``
  * ``result['gen'][i]['Qg']``

  hold the real (Pg) and reactive (Qg) power injection values at the optimal solution for generator `i`.
- If a feasible solution has been found, the additional ``branch`` entries

  * ``result["branch"][i]["Pf"]``
  * ``result["branch"][i]["Pt"]``
  * ``result["branch"][i]["Qf"]``
  * ``result["branch"][i]["Qt"]``
  * ``result["branch"][i]["switching"]``

  hold real (P) and reactive (Q) power injection values into the `from` (f) and into the `to` (t) end at the optimal solution point for branch `i`. The ``switching`` field holds the information whether a branch is turned on (1) or off (0) in the given result.


.. testcode:: opf

    result = opf.solve_opf(case, opftype="AC")

.. testoutput:: opf
    :hide:
    :options: +NORMALIZE_WHITESPACE +ELLIPSIS

    ...
    Optimize a model with 73 rows, 107 columns and 208 nonzeros
    ...

.. doctest:: opf
    :options: +NORMALIZE_WHITESPACE

    >>> result['success']
    1

    >>> result['bus'][0]
    {... 'Vm': 1.09..., 'Va': 0, ...}

    >>> result['branch'][1]
    {... 'Pf': 35.2..., 'Pt': -35.0..., 'Qf': -3.8..., 'Qt': -13.8..., ...}

    >>> result['gen'][2]
    {... 'Pg': 94.1..., 'Qg': -22.6..., ...}

We can see that the respective entries for ``bus`` and ``gen`` changed compared to the case dictionary, because they are different from the input at the optimal solution point. We also see that addtional fields have been created in the ``branch`` dictionary to hold solution information.


.. _recommended-label:

Recommended Literature
----------------------

Power systems and the optimal power flow problem are well studied. For a more comprehensive descrition, we recommend the following literature.


- G. Andersson. *Modelling and Analysis of Electric Power Systems*. Power Systems Laboratory,
  ETH Zürich, 2004.
- A.R. Bergen and V. Vittal. *Power Systems Analysis*. Prentice-Hall, 1999.
- D. Bienstock. *Electrical Transmission Systems Cascades and Vulnerability, an Operations Research
  viewpoint*. SIAM, 2015. ISBN 978-1-61197-415-7.
- D.K. Molzahn and I.A. Hiskens. *A survey of relaxations and approximations of the power flow
  equations*. Foundations and Trends in Electric Energy Systems, 4:1–221, 2019.
- J.D. Glover, M.S. Sarma, and T.J. Overbye. *Power System Analysis and Design*. CENGAGE
  Learning, 2012.
