import gurobipy as gp
from gurobipy import GRB


def set_gencost_objective(alldata, model):
    # Set quadratic objective using generator cost function coefficients

    GenPvar = alldata["LP"]["GenPvar"]
    gens = alldata["gens"]

    # Double checking: this should have been caught by input validation
    for gen in gens.values():
        assert gen.costdegree >= 1
        assert gen.costdegree == len(gen.costvector) - 1
        if gen.costdegree > 2:
            assert all(coeff == 0 for coeff in gen.costvector[:-3])

    # Quadratic terms are included directly (quadcostvar was never used)
    objective_quadratic = gp.quicksum(
        gen.costvector[-3] * GenPvar[gen] * GenPvar[gen]
        for gen in gens.values()
        if len(gen.costvector) >= 3
    )

    objective_linear = gp.quicksum(
        GenPvar[gen] * gen.costvector[-2] for gen in gens.values()
    )

    # The constant is added via a fixed variable
    constvar = model.addVar(lb=1.0, ub=1.0, name="constant")
    objective_constant = constvar * sum(gen.costvector[-1] for gen in gens.values())

    model.setObjective(
        objective_constant + objective_linear + objective_quadratic, sense=GRB.MINIMIZE
    )
