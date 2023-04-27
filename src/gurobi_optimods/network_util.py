import logging

import gurobipy as gp
import numpy as np
import scipy.sparse as sp
from gurobipy import GRB

try:
    import networkx as nx
except ImportError:
    nx = None

logger = logging.getLogger(__name__)


def _convert_pandas_to_digraph(edge_data, node_data):
    """
    Convert from a pandas DataFrame to a networkx.DiGraph with the appropriate
    attributes. For edges: `capacity`, and `cost`. For nodes: `demand`.
    """
    G = nx.from_pandas_edgelist(
        edge_data.reset_index(), create_using=nx.DiGraph(), edge_attr=True
    )
    for i, d in node_data.iterrows():
        G.add_node(i, demand=d.demand)
    return G


def _convert_pandas_to_scipy(edge_data, node_data):
    """
    Convert from a pandas DataFrame to several scipy.sparse.coo_matrix contain
    the graph structure, the capacity and cost values per edge, and the demand
    values per node.
    """
    coords = edge_data.index.to_numpy()

    a0 = np.array([c[0] for c in coords])
    a1 = np.array([c[1] for c in coords])

    data = np.ones(len(coords), dtype=np.int64)
    G = sp.coo_matrix((data, (a0, a1)))

    data = edge_data["capacity"].values
    cap = sp.coo_matrix((data, (a0, a1)))

    data = edge_data["cost"].values
    costs = sp.coo_matrix((data, (a0, a1)))

    return G, cap, costs, node_data["demand"].values


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
    indices = np.column_stack((edge_source, edge_target)).reshape(-1, order="C")
    indptr = np.arange(0, 2 * edge_source.shape[0] + 2, 2)
    ones = np.ones(edge_source.shape)
    data = np.column_stack((ones * -1.0, ones)).reshape(-1, order="C")
    A = sp.csc_array((data, indices, indptr))

    logger.info("Solving min-cost flow with {0} nodes and {1} edges".format(*A.shape))

    # Solve model with gurobi, return cost and flows
    with gp.Model(env=env) as m:
        m.ModelSense = GRB.MINIMIZE
        x = m.addMVar(A.shape[1], lb=0, ub=capacity, obj=cost)
        m.addMConstr(A, x, GRB.EQUAL, demand)
        m.optimize()
        if m.Status == GRB.INFEASIBLE:
            raise ValueError("Unsatisfiable flows")
        return m.ObjVal, x.X
