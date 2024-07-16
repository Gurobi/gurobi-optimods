import numpy as np


def _sort_key(x):
    return str(x)


def check_solution_pandas(solution, candidates):
    # Checks whether the solution (`pd.Series`) matches any of the list of
    # candidates (containing `dict`)
    if any(
        sorted(list(zip(solution.index.to_list(), solution.to_list())), key=_sort_key)
        == sorted(c, key=_sort_key)
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
    sol_list = sorted(
        [((i, j), data["flow"]) for i, j, data in solution.edges(data=True)],
        key=_sort_key,
    )
    if any(sol_list == sorted(c, key=_sort_key) for c in candidates):
        return True
    return False
