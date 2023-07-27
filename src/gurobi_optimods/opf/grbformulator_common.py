import gurobipy as gp
from gurobipy import GRB


def set_gencost_objective(alldata, model):
    # Set quadratic objective using generator cost function coefficients

    GenPvar = alldata["LP"]["GenPvar"]
    gens = alldata["gens"]

    # TODO fail earlier if not linear or quadratic.
    for gen in gens.values():
        assert gen.costdegree >= 2
        if gen.costdegree > 2:
            assert all(coeff == 0 for coeff in gen.costvector[:-3])

    # Quadratic parts are computed directly (quadcostvar was never used)
    # FIXME: wrong for NY case
    objective_quadratic = gp.quicksum(
        gen.costvector[0] * GenPvar[gen] * GenPvar[gen]
        for gen in gens.values()
        if gen.costdegree == 2 and gen.costvector[0] != 0
    )

    # Linear parts are added via an equality constraint
    lincostvar = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, name="lincost")
    objective_linear = gp.quicksum(
        GenPvar[gen] * gen.costvector[-2] for gen in gens.values()
    )
    model.addConstr(objective_linear == lincostvar, name="lincostdef")

    # The constant is added via a fixed variable
    constvar = model.addVar(lb=1.0, ub=1.0, name="constant")
    objective_constant = constvar * sum(gen.costvector[-1] for gen in gens.values())

    model.setObjective(
        objective_constant + lincostvar + objective_quadratic, sense=GRB.MINIMIZE
    )
