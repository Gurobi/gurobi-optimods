import numpy as np


def check_solution_pandas(solution, candidates):
    # Checks whether the solution (`pd.Series`) matches any of the list of
    # candidates (containing `pd.Series`)
    if any(solution.reset_index().equals(c.reset_index()) for c in candidates):
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
    for candidate in candidates:

        def edge_sort(row):
            return (str(row[0]), str(row[1]))

        if sorted(candidate, key=edge_sort) == sorted(
            list(solution.edges(data=True)), key=edge_sort
        ):
            return True
    return False
