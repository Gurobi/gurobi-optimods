Optimal Power Flow
==================

The operation of power systems relies on a number of optimization tasks, known as *optimal power flow* (*OPF*) problems. The objective of a standard OPF problem is to minimize operation cost such that the underlying grid constraints, such as generation, demand, and voltage limits, are satisfied.

In this mod, we consider the cases of *alternating current* (*AC*) and *direct current* (*DC*) OPF formulations. The ACOPF problem in its natural form requires the introduction of complex numbers in order to formulate the voltage. We circumvent the usage of complex numbers by using the so-called cartesian-coordinates formulation where we introduce additional optimization variables to formulate the complex terms via nonconvex quadratic relationships. The DCOPF problem is an approximation of the ACOPF problem where additional assumptions are made in order to make the optimization model linear. While the additional assumptions result in potential loss of solution accuracy, they make the DCOPF problem a lot easier to solve. This is especially useful if the solution accuracy can be neglected in favor of solution time and problem size. For additional information regarding the formulations, please refer to `More Complete Descriptions`_.

Here we assume basic familiarity with concepts such as *voltage* (potential energy), *current* (charge flow), and *power* (instantaneous energy generation or consumption).  The engineering community also uses the term *bus* (nodes in a network, with some simplification) and *branch* (arc in a network, a connection between two buses, typically a line or a transformer). For more details and more comprehensive descriptions of power systems and the underlying problems, please refer to the `Recommended Literature`_ section found further down below.


Problem Specification
---------------------

As input data for an OPF problem we have a power system (*power grid* for non-engineers) consisting of a network of buses and branches. Each branch is given as
an ordered pair :math:`km`, where :math:`k` is the *from* bus and :math:`m` is the *to* bus.

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

9. In modern versions discussed in the `Recommended Literature`_, network devices such as FACTS, phase-shifters, and impedance corrections are also incorporated into the model. We do not account for these in this mod.


More Complete Descriptions
--------------------------

.. toctree::
   :maxdepth: 1

   acopf
   dcopf


Solving an OPF Problem
----------------------

Input
~~~~~

This mod has multiple callable API and convenience functions. All of the main API functions take a so-called *case dictionary* as input. The case dictionary holds all essential information about the underlying network, i.e., information about buses, branch connections, and generators. The case dictionary is meant to follow the `MATPOWER Case Format conventions <https://matpower.org/docs/ref/matpower7.1/lib/caseformat.html>`_. In particular the case dictionary for this mod expects a dictionary with keys ``baseMVA``, ``bus``, ``branch``, ``gen``, and ``gencost``. All other entries of the dictionary are ignored in the current version of the mod. The value stored via the ``baseMVA`` key is a numerical float value. The values stored in the case dictionary via keys ``bus``, ``branch``, ``gen``, and ``gencost`` are lists of dictionaries, where each dictionary holds specific data about the particular object. Every single object is defined by a dictionary holding entries following the `MATPOWER Case Format <https://matpower.org/docs/ref/matpower7.1/lib/caseformat.html>`_, e.g., every bus has a bus number ``bus_i``, real power demand ``Pd``, etc. In the below example code, we read in a pre-defined case dictionary from our dataset.

.. doctest:: opf
    :options: +NORMALIZE_WHITESPACE

    >>> from gurobi_optimods.datasets import load_opfdictcase
    >>> case = load_opfdictcase()

    >>> case['baseMVA']
    100.0

    >>> buses = case['bus']
    >>> buses[1]
    {'bus_i': 1, 'type': 3, 'Pd': 0.0, 'Qd': 0.0, 'Gs': 0.0, 'Bs': 0.0, 'area': 0.0, 'Vm': 1.0, 'Va': 1.0, 'baseKV': 345.0, 'zone': 1.0, 'Vmax': 1.1, 'Vmin': 0.9}

    >>> branches = case['branch']
    >>> branches[2]
    {'fbus': 4, 'tbus': 5, 'r': 0.017, 'x': 0.092, 'b': 0.158, 'rateA': 250.0, 'rateB': 250.0, 'rateC': 250.0, 'ratio': 1.0, 'angle': 0.0, 'status': 1, 'angmin': -360.0, 'angmax': 360.0}

    >>> generators = case['gen']
    >>> generators[3]
    {'bus': 3, 'Pg': 85.0, 'Qg': 0.0, 'Qmax': 300.0, 'Qmin': -300.0, 'Vg': 1, 'mBase': 100, 'status': 1, 'Pmax': 270.0, 'Pmin': 10.0, 'Pc1': 0, 'Pc2': 0, 'Qc1min': 0, 'Qc1max': 0, 'Qc2min': 0, 'Qc2max': 0, 'ramp_agc': 0, 'ramp_10': 0, 'ramp_30': 0, 'ramp_q': 0, 'apf': 0}

    >>> generatorcosts = case['gencost']
    >>> generatorcosts[3]
    {'costtype': 2, 'startup': 3000, 'shutdown': 0, 'n': 3, 'costvector': [0.1225, 1, 335]}

There is also the convenience function :meth:`gurobi_optimods.opf.read_case_from_mat_file` which reads in a standard MATLAB ``.mat`` data file holding the network data. The data stored in the ``.mat`` file has to follow the `MATPOWER Case Format conventions <https://matpower.org/docs/ref/matpower7.1/lib/caseformat.html>`_ in order to be accepted by the function. Below, we generate a case dictionary from a ``.mat`` file containing network data for a small 9 bus network.

.. doctest:: opf
    :options: +NORMALIZE_WHITESPACE

    >>> from gurobi_optimods.datasets import load_caseopfmat
    >>> from gurobi_optimods.opf import read_case_from_mat_file
    >>> casefile = load_caseopfmat("9")
    >>> case = read_case_from_mat_file(casefile)

    >>> case['baseMVA']
    100

    >>> buses = case['bus']
    >>> buses[1]  # doctest:+ELLIPSIS
    {'bus_i': 1, ...}


Optimization Process
~~~~~~~~~~~~~~~~~~~~

After generating a case dictionary, we can solve an OPF problem defined by the given network data. For this task, we use the :meth:`gurobi_optimods.opf.solve_opf_model` function. We can define the type of the OPF problem that we want to solve by defining the ``opftype`` argument when calling the function. Currently, the available options are ``AC``, ``AC_relax``, and ``DC``.

- The ``AC`` setting solves an ACOPF problem defined by the given network data. The ACOPF problem is formulated as a nonconvex bilinear model as described in the :doc:`acopf`.

- The ``AC_relax`` setting solves a Second Order Cone (SOC) relaxation of the nonconvex bilinear ACOPF problem formulation defined by the given network data. The relaxation is constructed by dropping nonconvex bilinear terms but simultaneously keeping the convex JABR inequalities, see :doc:`acopf` for more details.

- The ``DC`` setting solves a DCOPF problem defined by the given network data. The DCOPF problem is a linear approximation of the ACOPF problem, see :doc:`dcopf` for more details.

The default value of the ``opftype`` argument is to solve an ``AC`` problem.

.. doctest:: opf
    :options: +NORMALIZE_WHITESPACE

    >>> from gurobi_optimods.datasets import load_opfdictcase
    >>> from gurobi_optimods.opf import solve_opf_model
    >>> case = load_opfdictcase()
    >>> result = solve_opf_model(case, opftype="AC")  # doctest:+ELLIPSIS
    All settings:
    ...
    Constructed ACOPF model with 107 variables and 73 constraints.
    ...
    Optimize a model with 73 rows, 107 columns and 208 nonzeros
    ...
    <BLANKLINE>
    Model Status: optimal.
    <BLANKLINE>
    Objective value = 5296....


Result
~~~~~~

We successfully solved an ACOPF problem and retrieved a so-called *result dictionary*. The result dictionary follows the same `MATPOWER Case Format conventions <https://matpower.org/docs/ref/matpower7.1/lib/caseformat.html>`_ as the case dictionary. However, in the result dictionary the some object entries are modified compared to the input case dictionary. These modified fields hold the solution values of the optimization. In some cases, there are also additional fields to store the solution information.

The following fields in the result dictionary are altered or added to store solution data.

- The field ``result['et']`` holds the runtime value of the whole solution process in seconds.
- The field ``result['success']`` defines whether at least one feasible solution has been found (1) or no feasible solution is available (0).
- The field ``result['f']`` holds the solution objective value (only valid if ``result['success'] == 1``).
- If a feasible solution has been found, the ``bus`` entries

  * ``result['bus'][i]['Vm']``
  * ``result['bus'][i]['Va']``

  store the voltage magnitude (Vm) and voltange angle (Va) values in the optimal solution for bus `i`.
- If we solved a DCOPF problem, additional fields ``result['bus'][i]['mu']`` holf the shadow prices for balance constraints at bus `i`.
- If a feasible solution has been found, the ``gen`` entries

  * ``result['gen'][i]['Pg']``
  * ``result['gen'][i]['Qg']``

  hold the real (Pg) and reactive (Qg) power injection values at the optimal solution for generator `i`.
- If a feasible solution has been found, the additional ``branch`` entries

  * ``result["branch"][i]["Pf"]``
  * ``result["branch"][i]["Pt"]``
  * ``result["branch"][i]["Qf"]``
  * ``result["branch"][i]["Qt"]``

  hold real (P) and reactive (Q) power injection values into the `from` (f) and into the `to` (t) end at the optimal solution point for branch `i`.

.. doctest:: opf
    :options: +NORMALIZE_WHITESPACE

    >>> from gurobi_optimods.datasets import load_opfdictcase
    >>> from gurobi_optimods.opf import solve_opf_model
    >>> case = load_opfdictcase()
    >>> result = solve_opf_model(case, opftype="AC")  # doctest:+SKIP
    >>> result['success']
    1

    >>> result['bus'][1]
    {... 'Vm': 1.09..., 'Va': 0, ...}

    >>> result['branch'][2]
    {... 'Pf': 35.2..., 'Pt': -35.0..., 'Qf': -3.8..., 'Qt': -13.8..., ...}

    >>> result['gen'][3]
    {... 'Pg': 94.1..., 'Qg': -22.6..., ...}

We can see that the respective entries for ``bus`` and ``gen`` changed compared to the case dictionary, because they are different from the input at the optimal solution point. We also see that addtional fields have been created in the ``branch`` dictionary to hold solution information.

It is possible to turn the result dictionary into a MATLAB ``.mat`` data file via the :meth:`turn_result_into_mat_file` function.

.. code-block:: python

    >>> turn_result_into_mat_file(result)


Branch Switching
----------------

TODO Waiting for a comprehensable example


Graphical Representation of Feasible Solutions
----------------------------------------------

In addition to solving an OPF problem, this mode also provides the possibility to plot the obtained result as a graphical representation of the network. Since there are already very involved graphical tools to represent OPF solutions, the graphical representation provided by this mod is very basic.  In order to use this functionality, it is necessary to install the ``plotly`` package.

Coordinate Information
~~~~~~~~~~~~~~~~~~~~~~

In order to plot a previously obtained result, an additional input of coordinates for all buses in the network is necessary. The coordinates have to be provided as a *coordinate dictionary*. The recommended way to generate a coordinate dictionary is to use a ``.csv`` holding all coordinate data. The ``.csv`` file holding the coordinate data has to follow the format

.. code-block::

   index(starting with 0), busID, busname, latitude, longitude
   0, 1, B1, 44.492, -73.208
   1, 2, B2, 41.271, -73.953
   ...

Once a ``.csv`` file holding bus coordinate information is available, we can use the :meth:`read_coords_from_csv_file` function to automatically generate a coordinate dictionary. In the following example we use the ``case9coords.csv`` file to generate a coordinate dictionary

.. doctest:: opf
    :options: +NORMALIZE_WHITESPACE

    >>> from gurobi_optimods.datasets import load_filepath
    >>> from gurobi_optimods.opf import read_coords_from_csv_file
    >>> coordsfile = load_filepath("case9coords.csv")
    >>> coords_dict = read_coords_from_csv_file(coordsfile)
    >>> coords_dict[1]
    (44.492, -73.208)

    >>> coords_dict[2]
    (41.271, -73.953)


Plotting the Result
~~~~~~~~~~~~~~~~~~~

After obtaining a result dictionary as discussed in `Solving an OPF Problem`_ and generating a coordinates dictionary, we can generate a :class:`plotly.graph_objects.Figure` object that can be displayed in, e.g., a browser window. Please note that it is required to install the ``plotly`` package to use this functionality. In the following, we solve the DCOPF problem for a given network data for the city of New York and plot the result.

.. code-block::

    >>> import plotly
    >>> from gurobi_optimods.datasets import load_caseNYopf, load_filepath
    >>> from gurobi_optimods.opf import solve_opf_model, read_case_from_mat_file, read_coords_from_csv_file
    >>> from gurobi_optimods.graphics import generate_opf_solution_figure
    >>> casefile = load_caseNYopf()
    >>> case = read_case_from_mat_file(casefile)
    >>> result = solve_opf_model(case, opftype='DC')

    >>> coordsfile = load_filepath("nybuses.csv")
    >>> coords_dict = read_coords_from_csv_file(coordsfile)
    >>> fig = generate_opf_solution_figure(case, coords_dict, solution)
    >>> fig.show()

.. image:: ../figures/opf.png

In the above image, you can see the power grid generated out of the given network data together with the coordinate information. The colored circles depict generators and the amount of power they generate

- Black bus: Power generation :math:`\leq 75` and load :math:`< 50`
- Blue bus: Power generation :math:`\leq 75` and load :math:`\geq 50`
- Purple bus: Power generation :math:`> 75`
- Orange bus: Power generation :math:`> 150`
- Red bus: Power generation :math:`> 500`


Violations for Pre-defined Voltage Values
-----------------------------------------

In practice, it is likely that we have voltage magnitudes and voltage angles for each bus at hand and would like to know whether these values are actually feasible within a given network. To tackle this, we can use the :meth:`compute_violations_from_given_voltages` function. This function takes a *voltage dictionary* and a case dictionary as arguments and is discussed in more detail below.


Voltage Information
~~~~~~~~~~~~~~~~~~~

In order to compute possible violations for given voltage data, an additional input of voltage information for all buses in the network is necessary. The voltage magnitudes (Vm) and voltage angles (Va) have to be provided as a *voltage dictionary*. The recommended way to generate a voltage dictionary is to use a ``.csv`` holding all voltage data. The ``.csv`` file holding the voltage data has to follow the format

.. code-block::

   index(starting with 0), busID, busname, Vm, Va
   0, 1, B1, 1.089026, 0.000000
   1, 2, B2, 1.099999, 20.552543
   ...

Once a ``.csv`` file holding voltage information for every bus is available, we can use the :meth:`read_voltages_from_csv_file` function to automatically generate a voltage dictionary. In the following example we use the ``case9volts.csv`` file to generate a voltage dictionary

.. doctest:: opf
    :options: +NORMALIZE_WHITESPACE

    >>> from gurobi_optimods.datasets import load_filepath
    >>> from gurobi_optimods.opf import read_voltages_from_csv_file
    >>> voltsfile = load_filepath("case9volts.csv")
    >>> volts_dict = read_voltages_from_csv_file(voltsfile)
    >>> volts_dict[1]
    (1.089026, 0.0)

    >>> volts_dict[2]
    (1.099999, 20.552543)


Checking for Violations
~~~~~~~~~~~~~~~~~~~~~~~

Once we have a voltage dictionary at hand, we can check for possible model violations by calling the :meth:`compute_violations_from_given_voltages` function. In addition to the verbose output, the function returns a *violations dictionary* which similar to the case dictionary follows the `MATPOWER Case Format <https://matpower.org/docs/ref/matpower7.1/lib/caseformat.html>`_. However, the violations dictionary has additional fields storing the violations for particular buses and branches.

The following fields in the violations dictionary are added to store violations data.

- ``violation['bus'][i]['Vmviol']`` Voltage magnitude violation at bus `i`
- ``violation['bus'][i]['Pviol']`` real power injection violation at bus `i`
- ``violation['bus'][i]['Qviol']`` reactive power injection violation at bus `i`
- ``violation['branch'][i]['limitviol']`` branch limit violation at branch `i`

.. doctest:: opf
    :options: +NORMALIZE_WHITESPACE

    >>> from gurobi_optimods.datasets import load_filepath, load_caseopfmat
    >>> from gurobi_optimods.opf import read_voltages_from_csv_file, read_case_from_mat_file, compute_violations_from_given_voltages
    >>> voltsfile = load_filepath("case9volts.csv")
    >>> volts_dict = read_voltages_from_csv_file(voltsfile)
    >>> casefile = load_caseopfmat("9")
    >>> case = read_case_from_mat_file(casefile)
    >>> violations = compute_violations_from_given_voltages(case, volts_dict)
    >>> violations['branch'][7]['limitviol']
    66.33435016796234

    >>> violations['bus'][4]['Pviol']
    -318.8997836192236

We can see that among others, the limit at branch 7 and the real power injection at bus 4 are violated.


Inspecting Violations Graphically
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Similar to generating a graphical representation of a feasible solution, it is also possible to generate a figure representing the violations within a given power grid. We can use the :meth:`generate_opf_violations_figure` to generate a :class:`plotly.graph_objects.Figure` object that can be displayed in, e.g., a browser window. Please note that it is required to install the ``plotly`` package to use this functionality. In addition to bus coordinates and case information, we also need to provide the violations dictionary to the function. In the following we use the previously discussed violated solution and plot the result.

.. code-block::

    >>> import plotly
    >>> from gurobi_optimods.datasets import load_filepath, load_caseopfmat
    >>> from gurobi_optimods.opf import read_coords_from_csv_file, read_voltages_from_csv_file, read_case_from_mat_file, compute_violations_from_given_voltages
    >>> from gurobi_optimods.graphics import generate_opf_solution_figure
    >>> voltsfile = load_filepath("case9volts.csv")
    >>> volts_dict = read_voltages_from_csv_file(voltsfile)
    >>> casefile = load_caseopfmat("9")
    >>> case = read_case_from_mat_file(casefile)
    >>> coordsfile = load_filepath("case9coords.csv")
    >>> coords_dict = read_coords_from_csv_file(coordsfile)
    >>> violations = compute_violations_from_given_voltages(case, volts_dict)
    >>> fig = generate_opf_violations_figure(case, coords, violations)
    >>> fig.show()

.. image:: ../figures/violations_opf.png

In the above image, you can see the power grid generated out of the given network data together with the coordinate and violation information. The red circles depict buses where the voltage magnitude or real or reactive power injections are violated. Red marked branches depict branches with violated limits.

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