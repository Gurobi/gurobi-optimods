import numpy as np


def check_solution_pandas(solution, candidates):
    # Checks whether the solution (`pd.Series`) matches any of the list of
    # candidates (containing `dict`)
    if any(
        solution.reset_index(drop=True).equals(c.reset_index(drop=True))
        for c in candidates
    ):
        return True
    return False


def check_solution_scipy(solution, candidates):
    # Checks whether the solution (`sp.sparray`) matches any of the list of
    # candidates (containing `np.ndarray`)
    arr = solution.toarray()
    if any(np.array_equal(arr, c) for c in candidates):
        return True
    return False


def check_solution_networkx(solution, candidates):
    # Checks whether the solution (`nx.DiGraph`) matches any of the list of
    # candidates (containing tuples dict `{(i, j): data}`)
    sol_dict = {(i, j): d for i, j, d in solution.edges(data=True)}
    if any(sol_dict == c for c in candidates):
        return True
    return False
