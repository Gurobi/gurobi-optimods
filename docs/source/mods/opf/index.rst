Optimal Power Flow
==================


Problem Specification
---------------------

The operation of power systems relies on a number of optimization tasks, known as 'optimal power flow' of OPF problems.  Here we assume basic familiarity with concepts such as 'voltage' (potential energy), 'current' (charge flow) and 'power' (instantaneous energy generation or consumption).  The engineering community also
uses the term 'bus' (nodes in a network, with some simplification) and 'branch' (a connection between two buses, typically a line or a transformer).

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


- :footcite:t:`bienstockpower`
      - Introduction to various power related problems with emphasis on modeling and optimization.

.. toctree::
   :maxdepth: 1

   acopf
   dcopf

Give examples of the various input data structures. These inputs should be fixed,
so use doctests where possible.

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
    from gurobi_optimods.datasets import load_opfsettings, load_caseopf


    # load path to settings file
    settingsfile = load_opfsettings()
    # read settings file and return a settings dictionary
    settings = read_settings_from_file(settingsfile)
    # load path to case file
    casefile = load_caseopf("9")
    # read case file and return a case dictionary
    case = read_case_from_file(casefile)
    # solve opf model and return a solution and the final objective value
    solution = solve_opf_model(settings, case)

..  A snippet of the Gurobi log output here won't show in the rendered page,
    but serves as a doctest to make sure the code example runs. The ... lines
    are meaningful here, they will match anything in the output test.

.. testoutput:: mod
    :hide:

    ...
    Optimize a model with 218 rows, 134 columns and 541 nonzeros
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

.. footbibliography::
