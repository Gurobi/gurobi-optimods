"""
Metro Map: Computing an Octilinear Graph Representation
-------------------------------------------------------
"""

import logging
import math

import gurobipy as gp
import pandas as pd

try:
    import networkx as nx
except ImportError:
    nx = None

try:
    import plotly.graph_objs as go
    from plotly.subplots import make_subplots
except ImportError:
    go = None

try:
    import plotly.io as pio

    pio.renderers.default = "browser"
except ImportError:
    pio = None


from gurobi_optimods.utils import optimod

logger = logging.getLogger(__name__)


@optimod()
def metromap(
    graph: nx.Graph,
    linepath_data: pd.DataFrame = None,
    include_planarity: bool = True,
    penalty_edge_directions=1,
    penalty_line_bends=1,
    penalty_distance=1,
    *,
    create_env,
):
    """Compute a metromap drawing, i.e., an octilinear network representation
     considering the origin geographical data if available given as geographic position of the nodes

    Parameters
    ----------
    graph : networkx graph with node attribute 'pos' being a tuple of x- and y-coordinate
    linepath_data : DataFrame
        DataFrame with information on the line routes/paths.
        It must include "linename", "edge_source", and "edge_target".
        This data frame could also be empty
    include_planarity : Bool
        parameter to turn off the consideration of planarity constraints
    penalty_edge_directions : Float
        weight for original direction part in the objective, default is 1
    penalty_line_bends : Float
        weight for line-bend part in the objective, default is 1
    penalty_distance : Float
        weight for distance part in the objective, default is 1

    Returns
    -------
    tuple
        - networkx graph with node property pos_oct representing the x and y coordinates of the octilinear representation
        - Python dict for direction for each edge (needed for plotting)
    """
    if nx is None:
        raise RuntimeError("Mod needs networkx package!")

    missing_data = False
    if linepath_data is not None:
        if "edge_source" not in linepath_data.columns:
            logger.info("column edge_source not present in linepath_data")
            missing_data = True
        if "edge_target" not in linepath_data.columns:
            logger.info("column edge_target not present in linepath_data")
            missing_data = True
        if linepath_data.isnull().values.any():
            logger.info("some value is nan or format is not correct for linepath_data")
            missing_data = True
    if missing_data:
        raise ValueError("Cannot run optimization. Some data is wrong or missing!")

    # octilinear graph representation needs a maximum node degree of 8
    # check node degree and abort if infeasible
    for v in graph.nodes:
        if graph.degree(v) > 8:
            logger.info(
                f"Node with number {v} has node degree {len(graph.edges(v))}. "
                f"Octilinear representation is not possible for node degree larger than 8"
            )
            return  ## value ???

    if linepath_data is None:
        linepaths = pd.DataFrame()
    else:
        linepaths = (
            linepath_data.set_index("linename")
            .groupby(["linename"])
            .apply(
                lambda x: [(k, v) for k, v in zip(x["edge_source"], x["edge_target"])]
            )
        )

    # check whether graph contains node position attribute
    pos = nx.get_node_attributes(graph, "pos")
    if len(pos) > 0:
        geodata = True
    else:
        geodata = False

    # remove isolated nodes
    graph.remove_nodes_from(list(nx.isolates(graph)))

    # restrict parameter to the interval of [0,100]
    if penalty_edge_directions < 0:
        penalty_edge_directions = 0
    if penalty_edge_directions > 100:
        penalty_edge_directions = 100
    if penalty_line_bends < 0:
        penalty_line_bends = 0
    if penalty_line_bends > 100:
        penalty_line_bends = 100
    if penalty_distance < 0:
        penalty_distance = 0
    if penalty_distance > 100:
        penalty_distance = 100

    return create_model(
        graph,
        linepaths,
        include_planarity,
        geodata,
        penalty_edge_directions,
        penalty_line_bends,
        penalty_distance,
        create_env,
    )


def create_model(
    graph,
    linepaths,
    include_planarity,
    geodata,
    penalty_edge_directions,
    penalty_line_bends,
    penalty_distance,
    create_env,
):
    """Create the model
    Parameters
    ----------
    graph : networkx graph with node attribute 'pos' being a tuple of x- and y-coordinate
    linepaths : DataFrame
        DataFrame with information on the line routes/paths.
        This data frame could also be empty
    include_planarity : Bool
        parameter to turn off the consideration of planarity constraints
    geofata : Bool
        if original positions were given
    penalty_edge_directions : Float
        weight for original direction part in the objective
    penalty_line_bends : Float
        weight for line-bend part in the objective
    penalty_distance : Float
        weight for distance part in the objective

    Returns
    -------
    tuple
        - networkx graph with node property pos_oct representing the x and y coordinates of the octilinear representation
        - Python dict for direction for each edge (needed for plotting)
    """

    edge_directions_out = dict()
    with create_env() as env, gp.Model(env=env) as model:
        num_nodes = len(graph.nodes)
        mindist = 1  # required distance in x- or y-coordinate for each edge
        # variables for position (coordinates)
        x = model.addVars(graph.nodes, ub=num_nodes, name="x-coord")
        y = model.addVars(graph.nodes, ub=num_nodes, name="y-coord")
        dist = {}
        for u, v in graph.edges:
            dist[u, v] = model.addVar(name=f"distance{u},{v}")
            dist[v, u] = model.addVar(name=f"distance{v},{u}")
            model.addConstr(dist[u, v] == dist[v, u])
            model.addConstr(dist[u, v] >= x[v] - x[u] - 1)
            model.addConstr(dist[u, v] >= y[v] - y[u] - 1)
            model.addConstr(dist[u, v] >= x[u] - x[v] - 1)
            model.addConstr(dist[u, v] >= y[u] - y[v] - 1)

        # direction for each edge (octilinear)
        edge_direction = {}
        for u, v in graph.edges:
            for i in range(8):
                edge_direction[u, v, i] = model.addVar(
                    vtype="B", name=f"direction[{u},{v},{i}]"
                )
                edge_direction[v, u, i] = model.addVar(
                    vtype="B", name=f"direction[{v},{u},{i}]"
                )
            # symmetric of back and forth direction (need to wait until all variables for (u,v) are created)
            model.addConstr(edge_direction[u, v, 0] == edge_direction[v, u, 4])
            model.addConstr(edge_direction[u, v, 1] == edge_direction[v, u, 5])
            model.addConstr(edge_direction[u, v, 2] == edge_direction[v, u, 6])
            model.addConstr(edge_direction[u, v, 3] == edge_direction[v, u, 7])
            model.addConstr(edge_direction[u, v, 4] == edge_direction[v, u, 0])
            model.addConstr(edge_direction[u, v, 5] == edge_direction[v, u, 1])
            model.addConstr(edge_direction[u, v, 6] == edge_direction[v, u, 2])
            model.addConstr(edge_direction[u, v, 7] == edge_direction[v, u, 3])

        # objective edge length
        obj = [0]
        for u, v in graph.edges:
            obj[0] += penalty_distance * dist[u, v]

        if geodata == True:
            # include original geographic information
            posOrig = nx.get_node_attributes(graph, "pos")
            # only allow neighboring positions
            create_orig_nodepostion_constraints(
                graph, edge_direction, obj, posOrig, penalty_edge_directions
            )
            # ensure the same edge order
            create_orig_edgeorder_constraints(graph, model, edge_direction, posOrig)

        # one direction for each edge
        model.addConstrs(
            (
                gp.quicksum(edge_direction[u, v, i] for i in range(8)) == 1
                for (u, v) in graph.edges
            ),
            name="OneDirection",
        )

        # restriction on x and y coordinates if a direction is chosen for an edge
        create_direction_constraints(graph, model, edge_direction, x, y, mindist)

        # define the bend for each two consecutive edges
        compute_bends(graph, model, edge_direction, obj, linepaths, penalty_line_bends)

        # add some helpful constraints (not needed for IP formulation)
        model.addConstrs(
            gp.quicksum(edge_direction[u, v, i] for (u, v) in graph.edges(w)) <= 1
            for w in graph.nodes
            if len(graph.edges(w)) > 1
            for i in range(8)
        )
        add_coordinateConstr(
            graph, model, edge_direction, x, y, mindist, num_nodes + mindist
        )

        model.setObjective(obj[0])
        model.update()

        if include_planarity:
            # helper variables for ensuring planarity
            startEdges = []
            gamma = {}
            for u1, v1 in graph.edges:
                startEdges.append((u1, v1))
                for u2, v2 in graph.edges:
                    if (u2, v2) in startEdges:
                        continue
                    if u1 == u2 or u1 == v2 or v1 == u2 or v1 == v2:
                        continue
                    for i in range(8):
                        gamma[u1, v1, u2, v2, i] = model.addVar(
                            vtype="b", name=f"planarity[{u1},{v1},{u2},{v2},{i}]"
                        )
                    objGamma = model.addVar(
                        lb=0, ub=1, name=f"objGamma[{u1},{v1},{u2},{v2}]"
                    )
                    model.addConstr(
                        gp.quicksum(gamma[u1, v1, u2, v2, i] for i in range(8))
                        + objGamma
                        == 1
                    )
                    obj[0] += 1000 * (objGamma)

            # necessary data for planarity callback
            model._mindist = mindist
            model._x = x
            model._y = y
            model._gamma = gamma
            model._bigM = num_nodes + num_nodes * mindist
            model._graph = graph
            model.setParam("LazyConstraints", 1)
            # optimize with callback for planarity constraints,
            # faster than adding all of them from the beginning
            model.optimize(planarity_callback)
        else:
            model.optimize()

        if model.solCount > 0:
            for v in graph.nodes:
                graph.add_node(v, pos_oct=(x[v].X, y[v].X))

            for u, v in graph.edges:
                for i in range(8):
                    if edge_direction[u, v, i].X > 0.5:
                        edge_directions_out[(u, v)] = i
                    if edge_direction[v, u, i].X > 0.5:
                        edge_directions_out[(v, u)] = i

    return graph, edge_directions_out


# Helper function
def counter_clockwise_sort(p, q):
    # Sort the positions given in q counter clockwise around p, start from east

    # Enumerate q to remember the original order
    q = list(enumerate(q))

    # Define a nested function for sorting
    def sort_by_angle(item):
        index, point = item
        dx, dy = point[0] - p[0], point[1] - p[1]
        angle = math.atan2(dy, dx)
        # Adjust the angle to start from East direction
        if angle < 0:
            angle += 2 * math.pi
        return angle

    # Sort the points counter-clockwise
    q.sort(key=sort_by_angle)

    # Return the sorted points along with their original indices
    return q


def create_orig_edgeorder_constraints(graph, model, edge_direction, posOrig):
    """Assume geographical node data is given.
        For each node in the graph order the neighbors counter-clock-wise and
        add constraints to ensure the preservation of this order in the octilinear
        graph representation

    Parameters
    ----------
    graph : graph (networkx)
        undirected networkx graph
    model : Gurobi model object
    edge_direction : tupledict (Gurobi variables)
    posOrig : node attribute (of graph)
    """
    for v in graph.nodes:
        # we can assume that the degree is at most 8
        if len(graph.edges(v)) <= 1:
            # nothing to do
            continue
        neighbors = []
        positions = []
        beta = model.addVars(len(graph.edges(v)), vtype="B", name=f"orientation{v}")
        model.addConstr(beta.sum() <= 1)
        for u, w in graph.edges(v):
            neighbors.append(w)
            positions.append(posOrig[w])
        sorted_q = counter_clockwise_sort(posOrig[v], positions)
        node = neighbors[sorted_q[0][0]]
        for i in range(1, len(sorted_q)):
            nextNode = neighbors[sorted_q[i][0]]
            model.addConstr(
                (beta[i - 1] == 0)
                >> (
                    gp.quicksum(k * edge_direction[v, node, k] for k in range(8))
                    <= gp.quicksum(k * edge_direction[v, nextNode, k] for k in range(8))
                    - 1
                )
            )

            node = nextNode
            if i == len(sorted_q) - 1:
                node0 = neighbors[sorted_q[0][0]]
                model.addConstr(
                    (beta[i] == 0)
                    >> (
                        gp.quicksum(k * edge_direction[u, node, k] for k in range(8))
                        <= gp.quicksum(
                            k * edge_direction[u, node0, k] for k in range(8)
                        )
                        - 1
                    )
                )


def create_orig_nodepostion_constraints(
    graph, edge_direction, obj, posOrig, penalty_edge_directions
):
    """Assume geographical node data is given.
        For each edge in the graph compute the direction from the original graph data.
        Fix variables such that for each edge only the original direction or the two adjacent
        directions are allowed.

    Parameters
    ----------
    graph : graph (networkx) - undirected networkx graph
    model : Gurobi model object
    edge_direction : tupledict (gurobi variables)
    obj : list - position 0 holds the objective funtion
    posOrig : node attribute (of graph) - original position
    penalty_edge_directions : penalty if not original position is chosen for the edge

    edge_direction : tupledict (gurobi variables)
    """
    for u, v in graph.edges:
        posU = posOrig[u]
        posV = posOrig[v]
        dx = posV[0] - posU[0]
        dy = posV[1] - posU[1]
        theta = math.atan2(dy, dx)
        theta_deg = math.degrees(theta)
        if theta_deg > -22.5 and theta_deg <= 22.5:
            # direction = "East"
            obj[0] += penalty_edge_directions * edge_direction[u, v, 1]
            edge_direction[u, v, 2].ub = 0
            edge_direction[u, v, 3].ub = 0
            edge_direction[u, v, 4].ub = 0
            edge_direction[u, v, 5].ub = 0
            edge_direction[u, v, 6].ub = 0
            obj[0] += penalty_edge_directions * edge_direction[u, v, 7]
        elif theta_deg > 22.5 and theta_deg <= 67.5:
            # direction = "North-East"
            obj[0] += penalty_edge_directions * edge_direction[u, v, 0]
            obj[0] += penalty_edge_directions * edge_direction[u, v, 2]
            edge_direction[u, v, 3].ub = 0
            edge_direction[u, v, 4].ub = 0
            edge_direction[u, v, 5].ub = 0
            edge_direction[u, v, 6].ub = 0
            edge_direction[u, v, 7].ub = 0
        elif theta_deg > 67.5 and theta_deg <= 112.5:
            # direction = "North"
            edge_direction[u, v, 0].ub = 0
            obj[0] += penalty_edge_directions * edge_direction[u, v, 1]
            obj[0] += penalty_edge_directions * edge_direction[u, v, 3]
            edge_direction[u, v, 4].ub = 0
            edge_direction[u, v, 5].ub = 0
            edge_direction[u, v, 6].ub = 0
            edge_direction[u, v, 7].ub = 0
        elif theta_deg > 112.5 and theta_deg <= 157.5:
            # direction = "North-West"
            edge_direction[u, v, 0].ub = 0
            edge_direction[u, v, 1].ub = 0
            obj[0] += penalty_edge_directions * edge_direction[u, v, 2]
            obj[0] += penalty_edge_directions * edge_direction[u, v, 4]
            edge_direction[u, v, 5].ub = 0
            edge_direction[u, v, 6].ub = 0
            edge_direction[u, v, 7].ub = 0
        elif theta_deg > 157.5 or theta_deg <= -157.5:
            # direction = "West"
            edge_direction[u, v, 0].ub = 0
            edge_direction[u, v, 1].ub = 0
            edge_direction[u, v, 2].ub = 0
            obj[0] += penalty_edge_directions * edge_direction[u, v, 3]
            obj[0] += penalty_edge_directions * edge_direction[u, v, 5]
            edge_direction[u, v, 6].ub = 0
            edge_direction[u, v, 7].ub = 0
        elif theta_deg > -157.5 and theta_deg <= -112.5:
            # direction = "South-West"
            edge_direction[u, v, 0].ub = 0
            edge_direction[u, v, 1].ub = 0
            edge_direction[u, v, 2].ub = 0
            edge_direction[u, v, 3].ub = 0
            obj[0] += penalty_edge_directions * edge_direction[u, v, 4]
            obj[0] += penalty_edge_directions * edge_direction[u, v, 6]
            edge_direction[u, v, 7].ub = 0
        elif theta_deg > -112.5 and theta_deg <= -67.5:
            # direction = "South"
            edge_direction[u, v, 0].ub = 0
            edge_direction[u, v, 1].ub = 0
            edge_direction[u, v, 2].ub = 0
            edge_direction[u, v, 3].ub = 0
            edge_direction[u, v, 4].ub = 0
            obj[0] += penalty_edge_directions * edge_direction[u, v, 5]
            obj[0] += penalty_edge_directions * edge_direction[u, v, 7]
        elif theta_deg > -67.5 and theta_deg <= -22.5:
            # direction = "South-East"
            obj[0] += penalty_edge_directions * edge_direction[u, v, 0]
            edge_direction[u, v, 1].ub = 0
            edge_direction[u, v, 2].ub = 0
            edge_direction[u, v, 3].ub = 0
            edge_direction[u, v, 4].ub = 0
            edge_direction[u, v, 5].ub = 0
            obj[0] += penalty_edge_directions * edge_direction[u, v, 6]


def create_direction_constraints(graph, model, edge_direction, x, y, mindist):
    """Define for each direction constraints on the position x and y for each edge

    Parameters
    ----------
    graph : graph (networkx) - undirected networkx graph
    model : Gurobi model object
    edge_direction : tupledict (gurobi variables)
    x : tupledict (gurobi variables) - x position of node
    y : tupledict (gurobi variables) - y position of node
    mindist : Scalar - minimum distance between two nodes
    """

    for u, v in graph.edges:
        model.addConstr((edge_direction[u, v, 0] == 1) >> (y[v] - y[u] == 0))
        model.addConstr((edge_direction[u, v, 0] == 1) >> (x[v] - x[u] >= mindist))
        model.addConstr((edge_direction[u, v, 1] == 1) >> (x[u] - y[u] == x[v] - y[v]))
        model.addConstr(
            (edge_direction[u, v, 1] == 1)
            >> (x[v] + y[v] - (x[u] + y[u]) >= 2 * mindist)
        )
        model.addConstr((edge_direction[u, v, 2] == 1) >> (x[u] == x[v]))
        model.addConstr((edge_direction[u, v, 2] == 1) >> (y[v] - y[u] >= mindist))
        model.addConstr((edge_direction[u, v, 3] == 1) >> (x[u] + y[u] == x[v] + y[v]))
        model.addConstr(
            (edge_direction[u, v, 3] == 1)
            >> (x[u] - y[u] - (x[v] - y[v]) >= 2 * mindist)
        )
        model.addConstr((edge_direction[u, v, 4] == 1) >> (y[u] == y[v]))
        model.addConstr((edge_direction[u, v, 4] == 1) >> (x[u] - x[v] >= mindist))
        model.addConstr((edge_direction[u, v, 5] == 1) >> (x[u] - y[u] == x[v] - y[v]))
        model.addConstr(
            (edge_direction[u, v, 5] == 1)
            >> (x[u] + y[u] - (x[v] + y[v]) >= 2 * mindist)
        )
        model.addConstr((edge_direction[u, v, 6] == 1) >> (x[u] == x[v]))
        model.addConstr((edge_direction[u, v, 6] == 1) >> (y[u] - y[v] >= mindist))
        model.addConstr((edge_direction[u, v, 7] == 1) >> (x[u] + y[u] == x[v] + y[v]))
        model.addConstr(
            (edge_direction[u, v, 7] == 1)
            >> (x[v] - y[v] - (x[u] - y[u]) >= 2 * mindist)
        )


def compute_bends(graph, model, edge_direction, obj, linepaths, penalty_line_bends):
    """Compute bends for each two consecutive lines. If there are linepaths covering
    both edges, add variables to model and objective that account for the "line bend".

    Parameters
    ----------
    graph : graph (networkx) - undirected networkx graph
    model : Gurobi model object
    edge_direction : tupledict (gurobi variables)
    obj : list - position 0 holds the objective funtion
    linepaths : pandas Dataframe - linepaths
    penalty_line_bends : penalty to account for bends of lines
    """
    bend = {}
    for v in graph.nodes:
        templist = list(graph.edges(v))
        for idx, (u1, w1) in enumerate(templist):
            # when connecting the first edge with all other edges, only one combination can have 0 bends
            linExpr0bend = 0
            for u2, w2 in templist[idx + 1 :]:
                # u1 = u2 and we have chain w1 - u1/2 - w2
                numlines = 0
                # compute number of lines covering both edges
                for i in linepaths.index:
                    if (w1, u1) in linepaths[i] and (u2, w2) in linepaths[i]:
                        numlines += 1
                    elif (w2, u2) in linepaths[i] and (u1, w1) in linepaths[i]:
                        numlines += 1
                if (
                    numlines > 0
                ):  # variables are only relevant when part of the objective
                    bend[w1, u1, w2, 0] = model.addVar(
                        vtype="B", name=f"bends[{w1},{u1},{w2},0]"
                    )
                    bend[w2, u1, w1, 0] = model.addVar(
                        vtype="B", name=f"bends[{w2},{u1},{w1},0]"
                    )
                    bend[w1, u1, w2, 1] = model.addVar(
                        vtype="B", name=f"bends[{w1},{u1},{w2},1]"
                    )
                    bend[w2, u1, w1, 1] = model.addVar(
                        vtype="B", name=f"bends[{w2},{u1},{w1},1]"
                    )
                    bend[w1, u1, w2, 2] = model.addVar(
                        vtype="B", name=f"bends[{w1},{u1},{w2},2]"
                    )
                    bend[w2, u1, w1, 2] = model.addVar(
                        vtype="B", name=f"bends[{w2},{u1},{w1},2]"
                    )
                    bend[w1, u1, w2, 3] = model.addVar(
                        vtype="B", name=f"bends[{w1},{u1},{w2},3]"
                    )
                    bend[w2, u1, w1, 3] = model.addVar(
                        vtype="B", name=f"bends[{w2},{u1},{w1},3]"
                    )
                    model.addConstr(
                        bend[w1, u1, w2, 0]
                        + bend[w1, u1, w2, 1]
                        + bend[w1, u1, w2, 2]
                        + bend[w1, u1, w2, 3]
                        == 1
                    )
                    model.addConstrs(
                        bend[w1, u1, w2, k] == bend[w2, u1, w1, k] for k in range(4)
                    )

                    model.addConstrs(
                        bend[w1, u1, w2, 0]
                        >= (edge_direction[w1, u1, k] + edge_direction[u2, w2, k] - 1)
                        for k in range(8)
                    )

                    model.addConstrs(
                        bend[w1, u1, w2, j]
                        >= (
                            edge_direction[w1, u1, k]
                            + edge_direction[u2, w2, k + j]
                            - 1
                        )
                        for j in range(1, 4)
                        for k in range(8 - j)
                    )
                    model.addConstrs(
                        bend[w1, u1, w2, j]
                        >= (
                            edge_direction[w1, u1, k]
                            + edge_direction[u2, w2, k - j]
                            - 1
                        )
                        for j in range(1, 4)
                        for k in range(j, 8)
                    )
                    model.addConstrs(
                        bend[w1, u1, w2, 8 - j]
                        >= (
                            edge_direction[w1, u1, k]
                            + edge_direction[u2, w2, k + j]
                            - 1
                        )
                        for j in range(5, 8)
                        for k in range(8 - j)
                    )
                    model.addConstrs(
                        bend[w1, u1, w2, 8 - j]
                        >= (
                            edge_direction[w1, u1, k]
                            + edge_direction[u2, w2, k - j]
                            - 1
                        )
                        for j in range(5, 8)
                        for k in range(j, 8)
                    )

                    obj[0] += (
                        penalty_line_bends
                        * numlines
                        * gp.quicksum(k * bend[w1, u1, w2, k] for k in range(4))
                    )
                    linExpr0bend += bend[w1, u1, w2, 0]

            if idx == 0:
                # this is a helpful cut but not necessary for the IP formulation
                model.addConstr(linExpr0bend <= 1, name=f"bend_{v}")


def planarity_callback(model, where):
    """Callback to handle planarity constraints"""
    if where == gp.GRB.Callback.MIPSOL:
        mindist = model._mindist
        mindistFeas = mindist - 1e-6
        bigM = model._bigM
        graph = model._graph
        x = model.cbGetSolution(model._x)
        y = model.cbGetSolution(model._y)
        startEdges = []
        for u1, v1 in graph.edges:
            startEdges.append((u1, v1))
            for u2, v2 in graph.edges:
                if (u2, v2) in startEdges:
                    continue
                if u1 == u2 or u1 == v2 or v1 == u2 or v1 == v2:
                    continue

                # check whether (u2, v2) is north of u1,v1 -> y coordinate increases
                if (
                    y[u2] - y[u1] >= mindistFeas
                    and y[u2] - y[v1] >= mindistFeas
                    and y[v2] - y[u1] >= mindistFeas
                    and y[v2] - y[v1] >= mindistFeas
                ):
                    continue
                # check if (u2, v2) is west of u1,v1 -> x coordinate decreases
                if (
                    x[u1] - x[u2] >= mindistFeas
                    and x[u1] - x[v2] >= mindistFeas
                    and x[v1] - x[u2] >= mindistFeas
                    and x[v1] - x[v2] >= mindistFeas
                ):
                    continue
                # check if (u2, v2) is south of u1,v1 -> y coordinate decreases
                if (
                    y[u1] - y[u2] >= mindistFeas
                    and y[u1] - y[v2] >= mindistFeas
                    and y[v1] - y[u2] >= mindistFeas
                    and y[v1] - y[v2] >= mindistFeas
                ):
                    continue
                # check if (u2, v2) is east of u1,v1 -> x coordinate inreases
                if (
                    x[u2] - x[u1] >= mindistFeas
                    and x[u2] - x[v1] >= mindistFeas
                    and x[v2] - x[u1] >= mindistFeas
                    and x[v2] - x[v1] >= mindistFeas
                ):
                    continue
                # check if (u2, v2) is north-west of u1,v1 -> z2 = (x-y) coordinate decreases
                if (
                    (x[u1] - y[u1]) - (x[u2] - y[u2]) >= mindistFeas
                    and (x[u1] - y[u1]) - (x[v2] - y[v2]) >= mindistFeas
                    and (x[v1] - y[v1]) - (x[u2] - y[u2]) >= mindistFeas
                    and (x[v1] - y[v1]) - (x[v2] - y[v2]) >= mindistFeas
                ):
                    continue
                # (u2, v2) is south-east of u1,v1 -> z2 = (x-y) coordinate inreases
                if (
                    (x[u2] - y[u2]) - (x[u1] - y[u1]) >= mindistFeas
                    and (x[u2] - y[u2]) - (x[v1] - y[v1]) >= mindistFeas
                    and (x[v2] - y[v2]) - (x[u1] - y[u1]) >= mindistFeas
                    and (x[v2] - y[v2]) - (x[v1] - y[v1]) >= mindistFeas
                ):
                    continue
                # (u2, v2) is south-east of u1,v1 -> z1 = (x+y) coordinate decreases
                if (
                    (x[u1] + y[u1]) - (x[u2] + y[u2]) >= mindistFeas
                    and (x[u1] + y[u1]) - (x[v2] + y[v2]) >= mindistFeas
                    and (x[v1] + y[v1]) - (x[u2] + y[u2]) >= mindistFeas
                    and (x[v1] + y[v1]) - (x[v2] + y[v2]) >= mindistFeas
                ):
                    continue
                # (u2, v2) is north-west of u1,v1 -> z1 = (x+y) coordinate inreases
                if (
                    (x[u2] + y[u2]) - (x[u1] + y[u1]) >= mindistFeas
                    and (x[u2] + y[u2]) - (x[v1] + y[v1]) >= -mindistFeas
                    and (x[v2] + y[v2]) - (x[u1] + y[u1]) >= mindistFeas
                    and (x[v2] + y[v2]) - (x[v1] + y[v1]) >= mindistFeas
                ):
                    continue

                # if we are here we found a violated planarity relation; add cut
                # (u2, v2) is north of u1,v1 -> y coordinate increases
                model.cbLazy(
                    model._y[u2] - model._y[u1]
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 2]) + mindist
                )
                model.cbLazy(
                    model._y[u2] - model._y[v1]
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 2]) + mindist
                )
                model.cbLazy(
                    model._y[v2] - model._y[u1]
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 2]) + mindist
                )
                model.cbLazy(
                    model._y[v2] - model._y[v1]
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 2]) + mindist
                )
                # (u2, v2) is west of u1,v1 -> x coordinate decreases
                model.cbLazy(
                    model._x[u1] - model._x[u2]
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 4]) + mindist
                )
                model.cbLazy(
                    model._x[u1] - model._x[v2]
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 4]) + mindist
                )
                model.cbLazy(
                    model._x[v1] - model._x[u2]
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 4]) + mindist
                )
                model.cbLazy(
                    model._x[v1] - model._x[v2]
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 4]) + mindist
                )
                # (u2, v2) is south of u1,v1 -> y coordinate decreases
                model.cbLazy(
                    model._y[u1] - model._y[u2]
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 6]) + mindist
                )
                model.cbLazy(
                    model._y[u1] - model._y[v2]
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 6]) + mindist
                )
                model.cbLazy(
                    model._y[v1] - model._y[u2]
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 6]) + mindist
                )
                model.cbLazy(
                    model._y[v1] - model._y[v2]
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 6]) + mindist
                )
                # (u2, v2) is east of u1,v1 -> x coordinate inreases
                model.cbLazy(
                    model._x[u2] - model._x[u1]
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 0]) + mindist
                )
                model.cbLazy(
                    model._x[u2] - model._x[v1]
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 0]) + mindist
                )
                model.cbLazy(
                    model._x[v2] - model._x[u1]
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 0]) + mindist
                )
                model.cbLazy(
                    model._x[v2] - model._x[v1]
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 0]) + mindist
                )
                # (u2, v2) is north-west of u1,v1 -> z2 = (x-y) coordinate decreases
                model.cbLazy(
                    (model._x[u1] - model._y[u1]) - (model._x[u2] - model._y[u2])
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 3]) + mindist
                )
                model.cbLazy(
                    (model._x[u1] - model._y[u1]) - (model._x[v2] - model._y[v2])
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 3]) + mindist
                )
                model.cbLazy(
                    (model._x[v1] - model._y[v1]) - (model._x[u2] - model._y[u2])
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 3]) + mindist
                )
                model.cbLazy(
                    (model._x[v1] - model._y[v1]) - (model._x[v2] - model._y[v2])
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 3]) + mindist
                )
                # (u2, v2) is south-east of u1,v1 -> z2 = (x-y) coordinate inreases
                model.cbLazy(
                    (model._x[u2] - model._y[u2]) - (model._x[u1] - model._y[u1])
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 7]) + mindist
                )
                model.cbLazy(
                    (model._x[u2] - model._y[u2]) - (model._x[v1] - model._y[v1])
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 7]) + mindist
                )
                model.cbLazy(
                    (model._x[v2] - model._y[v2]) - (model._x[u1] - model._y[u1])
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 7]) + mindist
                )
                model.cbLazy(
                    (model._x[v2] - model._y[v2]) - (model._x[v1] - model._y[v1])
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 7]) + mindist
                )
                # (u2, v2) is south-west of u1,v1 -> z1 = (x+y) coordinate decreases
                model.cbLazy(
                    (model._x[u1] + model._y[u1]) - (model._x[u2] + model._y[u2])
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 5]) + mindist
                )
                model.cbLazy(
                    (model._x[u1] + model._y[u1]) - (model._x[v2] + model._y[v2])
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 5]) + mindist
                )
                model.cbLazy(
                    (model._x[v1] + model._y[v1]) - (model._x[u2] + model._y[u2])
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 5]) + mindist
                )
                model.cbLazy(
                    (model._x[v1] + model._y[v1]) - (model._x[v2] + model._y[v2])
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 5]) + mindist
                )
                # (u2, v2) is north-east of u1,v1 -> z1 = (x+y) coordinate inreases
                model.cbLazy(
                    (model._x[u2] + model._y[u2]) - (model._x[u1] + model._y[u1])
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 1]) + mindist
                )
                model.cbLazy(
                    (model._x[u2] + model._y[u2]) - (model._x[v1] + model._y[v1])
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 1]) + mindist
                )
                model.cbLazy(
                    (model._x[v2] + model._y[v2]) - (model._x[u1] + model._y[u1])
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 1]) + mindist
                )
                model.cbLazy(
                    (model._x[v2] + model._y[v2]) - (model._x[v1] + model._y[v1])
                    >= -bigM * (1 - model._gamma[u1, v1, u2, v2, 1]) + mindist
                )


def add_coordinateConstr(graph, model, d, x, y, mindist, bigM):
    """Define additional cuts to improve the LP relaxation

    Parameters
    ----------
    graph : graph (networkx) - undirected networkx graph
    model : Gurobi model object
    edge_direction : tupledict (gurobi variables)
    x : tupledict (gurobi variables) - x position of node
    y : tupledict (gurobi variables) - y position of node
    mindist : Scalar - minimum distance between two nodes
    """

    for u, v in graph.edges:
        # one of the directions 0,1,2 then z1 increasing, y non-decreasing, x non-decreasing
        model.addConstr(
            x[v] + y[v]
            >= x[u]
            + y[u]
            + mindist
            - bigM * (d[u, v, 3] + d[u, v, 4] + d[u, v, 5] + d[u, v, 6] + d[u, v, 7])
        )
        model.addConstr(
            y[v]
            >= y[u]
            - bigM * (d[u, v, 3] + d[u, v, 4] + d[u, v, 5] + d[u, v, 6] + d[u, v, 7])
        )
        model.addConstr(
            x[v]
            >= x[u]
            - bigM * (d[u, v, 3] + d[u, v, 4] + d[u, v, 5] + d[u, v, 6] + d[u, v, 7])
        )
        # one of the directions 1,2,3 then y increasing, z2 non-increasing, z1 non-decreasing
        model.addConstr(
            y[v]
            >= y[u]
            + mindist
            - bigM * (d[u, v, 0] + d[u, v, 4] + d[u, v, 5] + d[u, v, 6] + d[u, v, 7])
        )
        model.addConstr(
            x[v] + y[v]
            >= x[u]
            + y[u]
            - bigM * (d[u, v, 0] + d[u, v, 4] + d[u, v, 5] + d[u, v, 6] + d[u, v, 7])
        )
        model.addConstr(
            x[v] - y[v]
            <= x[u]
            - y[u]
            + bigM * (d[u, v, 0] + d[u, v, 4] + d[u, v, 5] + d[u, v, 6] + d[u, v, 7])
        )
        # one of the directions 2,3,4 then z2 decreasing, y non-decreasing, x non-increasing
        model.addConstr(
            x[v] - y[v] + mindist
            <= x[u]
            - y[u]
            + bigM * (d[u, v, 0] + d[u, v, 1] + d[u, v, 5] + d[u, v, 6] + d[u, v, 7])
        )
        model.addConstr(
            y[v]
            >= y[u]
            - bigM * (d[u, v, 0] + d[u, v, 1] + d[u, v, 5] + d[u, v, 6] + d[u, v, 7])
        )
        model.addConstr(
            x[v]
            <= x[u]
            + bigM * (d[u, v, 0] + d[u, v, 1] + d[u, v, 5] + d[u, v, 6] + d[u, v, 7])
        )
        # one of the directions 3,4,5 then x decreasing, z1 non-increasing, z2 non-decreasing
        model.addConstr(
            x[v] + mindist
            <= x[u]
            + bigM * (d[u, v, 0] + d[u, v, 1] + d[u, v, 2] + d[u, v, 6] + d[u, v, 7])
        )
        model.addConstr(
            x[v] + y[v]
            <= x[u]
            + y[u]
            + bigM * (d[u, v, 0] + d[u, v, 1] + d[u, v, 2] + d[u, v, 6] + d[u, v, 7])
        )
        model.addConstr(
            x[v] - y[v]
            <= x[u]
            - y[u]
            + bigM * (d[u, v, 0] + d[u, v, 1] + d[u, v, 2] + d[u, v, 6] + d[u, v, 7])
        )
        # one of the directions 4,5,6 then z1 decreasing, y non-increasing, x non-decreasing
        model.addConstr(
            x[v] + y[v] + mindist
            <= x[u]
            + y[u]
            + bigM * (d[u, v, 0] + d[u, v, 1] + d[u, v, 2] + d[u, v, 3] + d[u, v, 7])
        )
        model.addConstr(
            y[v]
            <= y[u]
            + bigM * (d[u, v, 0] + d[u, v, 1] + d[u, v, 2] + d[u, v, 3] + d[u, v, 7])
        )
        model.addConstr(
            x[v]
            <= x[u]
            + bigM * (d[u, v, 0] + d[u, v, 1] + d[u, v, 2] + d[u, v, 3] + d[u, v, 7])
        )
        # one of the directions 5,6,7 then y decreasing, z2 non-decreasing, z1 non-decreasing
        model.addConstr(
            y[v] + mindist
            <= y[u]
            + bigM * (d[u, v, 0] + d[u, v, 1] + d[u, v, 2] + d[u, v, 3] + d[u, v, 4])
        )
        model.addConstr(
            x[v] - y[v]
            >= x[u]
            - y[u]
            - bigM * (d[u, v, 0] + d[u, v, 1] + d[u, v, 2] + d[u, v, 3] + d[u, v, 4])
        )
        model.addConstr(
            x[v] + y[v]
            <= x[u]
            + y[u]
            + bigM * (d[u, v, 0] + d[u, v, 1] + d[u, v, 2] + d[u, v, 3] + d[u, v, 4])
        )
        # one of the directions 6,7,0 then z2 increasing, x non-decreasing, y non-increasing
        model.addConstr(
            x[v] - y[v]
            >= x[u]
            - y[u]
            + mindist
            - bigM * (d[u, v, 1] + d[u, v, 2] + d[u, v, 3] + d[u, v, 4] + d[u, v, 5])
        )
        model.addConstr(
            x[v]
            >= x[u]
            - bigM * (d[u, v, 1] + d[u, v, 2] + d[u, v, 3] + d[u, v, 4] + d[u, v, 5])
        )
        model.addConstr(
            y[v]
            <= y[u]
            + bigM * (d[u, v, 1] + d[u, v, 2] + d[u, v, 3] + d[u, v, 4] + d[u, v, 5])
        )
        # one of the directions 7,0,1 then x increasing, z2 non-decreasing, z1 non-decreasing
        model.addConstr(
            x[v]
            >= x[u]
            + mindist
            - bigM * (d[u, v, 2] + d[u, v, 3] + d[u, v, 4] + d[u, v, 5] + d[u, v, 6])
        )
        model.addConstr(
            x[v] - y[v]
            >= x[u]
            - y[u]
            - bigM * (d[u, v, 2] + d[u, v, 3] + d[u, v, 4] + d[u, v, 5] + d[u, v, 6])
        )
        model.addConstr(
            x[v] + y[v]
            >= x[u]
            + y[u]
            - bigM * (d[u, v, 2] + d[u, v, 3] + d[u, v, 4] + d[u, v, 5] + d[u, v, 6])
        )


def plot_network(graph: nx.graph, directions: dict(), linepath_data: pd.DataFrame):
    """Plotting function for a metro map, the lines can be turned off and on via clicking in the legend

    Parameters
    ----------
    graph : networkx graph with node attribute 'pos_oct' being a tuple of x- and y-coordinate
    directions: dictionary providing the assigned direction (0 to 7) for each edge
    linepath_data : DataFrame
        DataFrame with information on the line routes/paths.
    Returns
    -------
    A browser opens with the plot
    """
    if pio is None or go is None:
        raise RuntimeError(
            "Plot not possible: plotly.graph_objs, plotly.subplots, and plotly.io are required for plotting the map"
        )

    if graph.number_of_nodes() == 0:
        logger.info(f"No nodes in graph, nothing to draw.")
        return
    pos_oct = nx.get_node_attributes(graph, "pos_oct")
    if len(pos_oct) == 0:
        logger.info(f"No node attribute with name pos_oct is given; cannot draw graph.")
        return

    if linepath_data is None or len(directions) == 0:
        logger.info(f"No lines given or no edge directions; nothing to draw.")
        return

    try:
        linepaths = (
            linepath_data.set_index("linename")
            .groupby(["linename"])
            .apply(
                lambda x: [(k, v) for k, v in zip(x["edge_source"], x["edge_target"])]
            )
        )
    except:
        raise ValueError(f"linepath_data does not contain relevant information")
    plot_lines(graph, linepaths, directions)


# plotting function #
def plot_lines(graph, linepaths, directions):
    num_nodes = graph.number_of_nodes()
    scale = 0.01
    if num_nodes > 30:
        scale = 0.02
    if num_nodes > 70:
        scale = 0.03
    if num_nodes > 110:
        scale = 0.04
    if num_nodes > 150:
        scale = 0.05

    # Create a subplot
    fig = make_subplots(rows=1, cols=1)

    # shift lines if they use the same edge
    edge_shift = dict()
    # initialize with empty list
    for u, v in graph.edges:
        edge_shift[(u, v)] = list()
        edge_shift[(v, u)] = edge_shift[(u, v)]

    for l in linepaths.index:
        # combine all edges of a line that have the same direction, then shift the whole bunch
        last_direction = -1
        combined_edges = []
        # only draw line for each edge once (in case the data contain back and forth direction of edges)
        covered_edges = []
        xcoord = []
        ycoord = []
        for ut, vt in linepaths[l]:
            # we detect a an edge that was considered before for this line
            if (vt, ut) in covered_edges or (ut, vt) in covered_edges:
                # handle edges of the last bunch
                if len(combined_edges) > 0:
                    add_coordinates(
                        graph,
                        combined_edges,
                        edge_shift,
                        directions,
                        xcoord,
                        ycoord,
                        scale,
                    )
                # reset information for edge-bunch
                last_direction = -1
                combined_edges = []
                continue
            covered_edges.append((ut, vt))
            if last_direction == -1:
                last_direction = directions[(ut, vt)]
                combined_edges.append((ut, vt))
            else:
                if last_direction == directions[(ut, vt)]:
                    combined_edges.append((ut, vt))
                else:  # new direction, handle line for combined edges (=last bunch)
                    add_coordinates(
                        graph,
                        combined_edges,
                        edge_shift,
                        directions,
                        xcoord,
                        ycoord,
                        scale,
                    )
                    last_direction = directions[(ut, vt)]
                    combined_edges = [(ut, vt)]
        # handle last bunch of edges
        if len(combined_edges) > 0:
            add_coordinates(
                graph, combined_edges, edge_shift, directions, xcoord, ycoord, scale
            )
        # draw all edges for the line, make the whole line selectable
        fig.add_trace(
            go.Scatter(x=xcoord, y=ycoord, mode="lines", line_width=3, name=l)
        )

    # plot all nodes, nodes are not selectable
    for n in graph.nodes:
        (x, y) = graph.nodes[n]["pos_oct"]
        fig.add_trace(
            go.Scatter(
                x=[x],
                y=[y],
                mode="markers",
                marker=dict(
                    color="white",
                    size=math.ceil(1 / math.sqrt(num_nodes) * 100),
                    line=dict(width=2, color="black"),
                ),
                showlegend=False,
            )
        )

    # do not show axis, prefer white background
    fig.update_layout(xaxis_visible=False, yaxis_visible=False)
    fig.update_layout(plot_bgcolor="white")

    # Show the figure
    fig.show()


# helper function for plotting #
def check_shift(combined_edges, edge_shift, shift_val):
    # check if this shift value was already selected for one of the edges
    for c in combined_edges:
        if shift_val in edge_shift[c]:
            return True
    return False


# helper function for plotting #
def shift_combined_edges(graph, combined_edges, edge_shift):
    # find the right shift, do not more than 10 shifts
    for i in [0, -1, 1, -2, 2, -3, 3, -4, 4, -5, 5]:
        check = check_shift(combined_edges, edge_shift, i)
        if check == False:
            return i
    return i


# helper function for plotting #
def add_coordinates(
    graph, combined_edges, edge_shift, directions, xcoord, ycoord, scale
):
    # compute x, and y coordinates for the edges in combined_edges
    # shift them if necessary
    shift = shift_combined_edges(graph, combined_edges, edge_shift)
    for u, v in combined_edges:
        (x1, y1) = graph.nodes[u]["pos_oct"]
        (x2, y2) = graph.nodes[v]["pos_oct"]
        cnt = shift * scale
        if directions[(u, v)] == 0 or directions[(u, v)] == 4:
            y1 += 1.5 * cnt
            y2 += 1.5 * cnt
        elif directions[(u, v)] == 1 or directions[(u, v)] == 5:
            x1 -= cnt
            x2 -= cnt
            y1 += cnt
            y2 += cnt
        elif directions[(u, v)] == 2 or directions[(u, v)] == 6:
            x1 -= 1.5 * cnt
            x2 -= 1.5 * cnt
        elif directions[(u, v)] == 3 or directions[(u, v)] == 7:
            x1 += cnt
            x2 += cnt
            y1 += cnt
            y2 += cnt
        xcoord.append(x1)
        xcoord.append(x2)
        ycoord.append(y1)
        ycoord.append(y2)
        edge_shift[(u, v)].append(shift)
