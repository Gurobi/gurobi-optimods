def load_case9branchswitching():
    # we alter the original case dictionary in order to
    # create an artifical case, where turning off 2 branches
    # produces a better solution
    from gurobi_optimods.datasets import load_opf_example

    casefile_dict = load_opf_example("case9")
    casefile_dict["branch"].extend(
        [
            {
                "fbus": 1,
                "tbus": 2,
                "r": 0.0,
                "x": 0.0576,
                "b": 0.0,
                "rateA": 250.0,
                "rateB": 250.0,
                "rateC": 250.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "angmin": -360.0,
                "angmax": 360.0,
            },
            {
                "fbus": 1,
                "tbus": 3,
                "r": 0.0,
                "x": 0.0576,
                "b": 0.0,
                "rateA": 250.0,
                "rateB": 250.0,
                "rateC": 250.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "angmin": -360.0,
                "angmax": 360.0,
            },
            {
                "fbus": 2,
                "tbus": 3,
                "r": 0.0,
                "x": 0.0576,
                "b": 0.0,
                "rateA": 250.0,
                "rateB": 250.0,
                "rateC": 250.0,
                "ratio": 1.0,
                "angle": 0.0,
                "status": 1,
                "angmin": -360.0,
                "angmax": 360.0,
            },
        ]
    )
    casefile_dict["gencost"][2]["costvector"] = [0.85, 10.2, 1200]

    return casefile_dict
