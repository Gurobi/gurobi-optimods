import gurobipy as gp
import numpy as np
import scipy.sparse as sp
from gurobipy import GRB


def solve_min_cost_flow(
    env: gp.Env,
    edge_source: np.ndarray,
    edge_target: np.ndarray,
    capacity: np.ndarray,
    cost: np.ndarray,
    demand: np.ndarray,
):
    """Solve a min-cost flow model for the given network.

    Formulated as min-cost flow (A, b, c, u) using the matrix-friendly API, where
    A = incidence matrix, b = demand, c = cost, u = capacity::

        min   c^T x
        s.t.  A x = b
              0 <= x <= u

    :param env: Gurobi environment to use for model
    :type env: :class:`gurobipy.Env`
    :param edge_source: Index array; ``edge_source[i]`` is the source node
        :math:`j` of edge :math:`i`
    :type edge_source: :class:`numpy.ndarray`
    :param edge_target: Index array; ``edge_target[i]`` is the source node
        :math:`j` of edge :math:`i`
    :type edge_target: :class:`numpy.ndarray`
    :param capacity: Capacity array; ``capacity[i]`` is the capacity of
        edge :math:`i`
    :type capacity: :class:`numpy.ndarray`
    :param cost: Cost array; ``cost[i]`` is the cost of edge :math:`i` (to
        be minimized)
    :type cost: :class:`numpy.ndarray`
    :param demand: Demand array; ``demand[j]`` is the net demand of node
        :math:`j`
    :type demand: :class:`numpy.ndarray`
    """

    # Create incidence matrix from edge lists.
    # FIXME: feels like a long way round to CSR?
    s = edge_source.shape[0]
    i = np.concatenate([edge_source, edge_target])
    j = np.concatenate([np.arange(s), np.arange(s)])
    data = np.ones(s * 2, dtype=float)
    data[:s] *= -1.0
    A = sp.coo_array((data, (i, j))).tocsr()

    # Balance, capacities and costs are given, just check shape
    # FIXME: throw informative errors
    nodes, edges = A.shape
    assert demand.shape == (nodes,)
    assert capacity.shape == (edges,)
    assert cost.shape == (edges,)

    # Add names if exporting for debugging
    dump_model = False
    if dump_model:
        names = np.array([f"x[{i},{j}]" for i, j in zip(edge_source, edge_target)])
    else:
        names = None

    # Solve model with gurobi, return cost and flows
    with gp.Model(env=env) as m:
        m.ModelSense = GRB.MINIMIZE
        x = m.addMVar(edges, lb=0, ub=capacity, obj=cost, name=names)
        m.addMConstr(A, x, GRB.EQUAL, demand)
        if dump_model:
            m.write("netflow.lp")
        m.optimize()
        if m.Status == GRB.INFEASIBLE:
            raise ValueError("Unsatisfiable flows")
        if dump_model:
            m.write("netflow.sol")
        return m.ObjVal, x.X
