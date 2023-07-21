"""

- TODO would prefer `datasets` had a function directly returing the coordinate
  and voltage data needed for examples (rather than getting a file path then
  calling the csv reader) but I can't think what to call it yet. For case9, the
  content is short enough that we could probably put the dictionary directly in
  the docs.

Usage examples:


    # Load and solve 9 bus example

    from gurobi_optimods import opf
    from gurobi_optimods.datasets import load_opf_example, load_opf_extra

    case = load_opf_example("case9")
    solution = opf.solve_opf(case, opftype="ac")


    # Plot results with given coordinates

    coordinates = load_opf_extra("case9-coordinates")
    figure = opf.solution_plot(case, coordinates, solution)
    figure.show()


    # Check violations and plot

    voltages = load_opf_extra("case9-voltages")
    violations = opf.compute_violations(case, voltages=voltages)
    figure = opf.violation_plot(case, coordinates, violations)
    figure.show()


    # Load a case from a MATPOWER .mat file

    case = opf.read_case_matpower("mycase.mat")


    # Load the New York example

    case = load_opf_example("caseNY")


"""

from gurobi_optimods.opf.api import compute_violations  # noqa: F401
from gurobi_optimods.opf.api import solve_opf  # noqa: F401
from gurobi_optimods.opf.graphics import solution_plot  # noqa: F401
from gurobi_optimods.opf.graphics import violation_plot  # noqa: F401
from gurobi_optimods.opf.io import read_case_matpower  # noqa: F401
from gurobi_optimods.opf.io import write_case_matpower  # noqa: F401
