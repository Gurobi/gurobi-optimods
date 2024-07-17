import numpy as np


def check_solution_pandas(solution, candidates):
    # Checks whether the solution (`pd.Series`) matches any of the list of
    # candidates. Each candidate is a list of tuples ((i, j), v) tuples,
    # compare with the solution in sorted order.
    solution_list = sorted(solution.items())
    return any(solution_list == sorted(candidate) for candidate in candidates)


def check_solution_scipy(solution, candidates):
    # Checks whether the solution (`sp.sparray`) matches any of the list of
    # candidates (containing `np.ndarray`)
    arr = solution.toarray()
    return any(np.array_equal(arr, c) for c in candidates)


def check_solution_networkx(solution, candidates):
    # Checks whether the solution (`nx.DiGraph`) matches any of the list of
    # candidates (containing tuples dict `{(i, j): data}`)
    solution_list = sorted(
        [((i, j), data["flow"]) for i, j, data in solution.edges(data=True)],
    )
    return any(solution_list == sorted(candidate) for candidate in candidates)
