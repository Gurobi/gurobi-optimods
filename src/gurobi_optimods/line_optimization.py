"""
Line Optimization in Public Transportation Networks
---------------------------------------------------
"""

import logging

import gurobipy as gp
import pandas as pd

try:
    import networkx as nx
except ImportError:
    nx = None

try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
except ImportError:
    mpl = None


from gurobi_optimods.utils import optimod

logger = logging.getLogger(__name__)


@optimod()
def line_optimization(
    node_data: pd.DataFrame,
    edge_data: pd.DataFrame,
    line_data: pd.DataFrame,
    linepath_data: pd.DataFrame,
    demand_data: pd.DataFrame,
    frequencies: list,
    shortest_paths: bool = True,
    *,
    create_env,
):
    """Solve the line planning problem.

    Parameters
    ----------
    node_data : DataFrame
        DataFrame with information on the nodes/stations. The frame must include "source".
    edge_data : DataFrame
        DataFrame with edges / connections between stations and associated attributes.
        The frame must include "source", "target", and "time"
    demand_data : DataFrame
        DataFrame with node/station demand information.
        It must include "source", "target", and "demand". The demand value must be non-negative.
    line_data : DataFrame
        DataFrame with general line information.
        It must include "linename", "capacity", "fix_cost", and "operating_cost".
    linepath_data : DataFrame
        DataFrame with information on the line routes/paths.
        It must include "linename", "edge_source", and "edge_target".
    frequency: List
        List with possible frequencies: How often the line can be operated in the considered
        time horizon.

    Returns
    -------
    tuple
        - Cost of the optimal line concept (line frequency)
        - list of line-frequency tuples (optimal line concept)
    """
    # check for missing or wrong data
    missing_data = False
    if "number" not in node_data.columns:
        logger.info("column nr not present in node_data")
        missing_data = True
    if "time" not in edge_data.columns:
        logger.info("column time not present in edge_data")
        missing_data = True
    if "source" not in edge_data.columns:
        logger.info("column source not present in edge_data")
        missing_data = True
    if "target" not in edge_data.columns:
        logger.info("column target not present in edge_data")
        missing_data = True
    if "linename" not in line_data.columns:
        logger.info("column linename not present in line_data")
        missing_data = True
    if "capacity" not in line_data.columns:
        logger.info("column capacity not present in line_data")
        missing_data = True
    if "fix_cost" not in line_data.columns:
        logger.info("column fix_cost not present in line_data")
        missing_data = True
    if "operating_cost" not in line_data.columns:
        logger.info("column operating_cost not present in line_data")
        missing_data = True
    if "linename" not in linepath_data.columns:
        logger.info("column linename not present in linepath_data")
        missing_data = True
    if "edge_source" not in linepath_data.columns:
        logger.info("column edge_source not present in linepath_data")
        missing_data = True
    if "edge_target" not in linepath_data.columns:
        logger.info("column edge_target not present in linepath_data")
        missing_data = True
    if "source" not in demand_data.columns:
        logger.info("column source not present in demand_data")
        missing_data = True
    if "target" not in demand_data.columns:
        logger.info("column target not present in demand_data")
        missing_data = True
    if "demand" not in demand_data.columns:
        logger.info("column demand not present in demand_data")
        missing_data = True
    if node_data.isnull().values.any():
        logger.info("some value is nan or format is not correct for node_data")
        missing_data = True
    if edge_data.isnull().values.any():
        logger.info("some value is nan or format is not correct for edge_data")
        missing_data = True
    if line_data.isnull().values.any():
        logger.info("some value is nan or format is not correct for line_data")
        missing_data = True
    if linepath_data.isnull().values.any():
        logger.info("some value is nan or format is not correct for linepath_data")
        missing_data = True
    if demand_data.isnull().values.any():
        logger.info("some value is nan or format is not correct for demand_data")
        missing_data = True
    if (demand_data < 0).values.any():
        logger.info(
            "All demand should be non-negative! Some value is negative in demand_data"
        )
        missing_data = True
    if missing_data:
        raise ValueError("Cannot run optimization. Some data is wrong or missing!")

    # prepare data
    nodes = node_data.set_index("number").to_dict("index")
    edges = edge_data.set_index(["source", "target"]).to_dict("index")
    lines = line_data.set_index("linename").to_dict("index")
    linepaths = (
        linepath_data.set_index("linename")
        .groupby(["linename"])
        .apply(lambda x: [(k, v) for k, v in zip(x["edge_source"], x["edge_target"])])
    )
    demands = demand_data.set_index(["source", "target"]).to_dict()["demand"]

    # remove edges that are not covered by any lines
    for k in list(edges.keys()):
        found = False
        for l in linepaths.index:
            if k in linepaths[l]:
                found = True
                break
        if found == False:
            del edges[k]

    if shortest_paths:
        if nx is None:
            logger.info(
                "Networkx is needed for strategy 1 but not available. Using strategy 2 instead."
            )
            return allow_all_paths(
                nodes,
                edges,
                lines,
                linepaths,
                demands,
                demand_data,
                frequencies,
                create_env,
            )
        return all_shortest_paths(
            nodes, edges, edge_data, lines, linepaths, demands, frequencies, create_env
        )
    else:
        return allow_all_paths(
            nodes,
            edges,
            lines,
            linepaths,
            demands,
            demand_data,
            frequencies,
            create_env,
        )


def all_shortest_paths(
    nodes, edges, edge_data, lines, linepaths, demands, frequencies, create_env
):
    """
    Strategy 1:
    - compute all shortest path for each OD pair using networkx all_shortest_paths algorithm
    - create passenger variable for each path
    """
    logger.info(
        "Starting line optimization using strategy 1 - allow only shortest paths."
    )
    G = nx.from_pandas_edgelist(
        edge_data, create_using=nx.DiGraph(), edge_attr=["time"]
    )

    with create_env() as env, gp.Model(env=env) as model:
        # add variables for lines and frequencies
        x = {}
        numlines = 0
        for l in linepaths.index:
            numlines += 1
            for f in frequencies:
                x[l, f] = model.addVar(
                    vtype=gp.GRB.BINARY,
                    obj=f * lines[l]["operating_cost"] + lines[l]["fix_cost"],
                    name=str(l) + str(f),
                )

        logger.info(
            f"Consider network with {len(nodes)} nodes, "
            f"{len(edges)} edges, and "
            f"{numlines} potential line routes"
        )

        # prepare capacity constraints
        cap = model.addConstrs(
            (
                (
                    gp.quicksum(
                        x[l, f] * lines[l]["capacity"] * f
                        for l in linepaths.index
                        for f in frequencies
                        if (u, v) in linepaths[l]
                    )
                )
                >= 0
                for (u, v) in edges
            ),
            name="capacity",
        )

        # compute all shortest paths (same travel time) for each OD pair
        # define a variable for each such shortest path
        f = {}  # stuv from s to t using edge uv
        for s, t in demands:
            if s not in G:
                raise ValueError(f"demand node {s} not found in edges")
            if t not in G:
                raise ValueError(f"demand node {t} not found in edges")
            try:
                paths = list(
                    nx.all_shortest_paths(G, source=s, target=t, weight="time")
                )
            except nx.NetworkXNoPath:
                raise ValueError(f"no path found for connection from {s} to {t}")
            demand_expr = 0
            for p in paths:
                y = model.addVar(
                    ub=demands[s, t],
                    vtype=gp.GRB.CONTINUOUS,
                    name="passpath" + str(s) + "," + str(t) + str(p),
                )
                demand_expr += y
                for node in range(len(p) - 1):
                    model.chgCoeff(cap[p[node], p[node + 1]], y, -1.0)
            model.addConstr(demand_expr == demands[s, t])

        # each line with at most one frequency
        model.addConstrs(
            (gp.quicksum(x[i, j] for j in frequencies) <= 1 for i in linepaths.index),
            name="one-freq",
        )

        model.optimize()

        # prepare return values
        obj_cost = -1
        lines_out = []
        if model.Status == gp.GRB.OPTIMAL:
            obj_cost = model.objVal
            for i in linepaths.index:
                for j in frequencies:
                    if x[i, j].X > 0.5:
                        lines_out.append((i, j))

        return obj_cost, lines_out


def allow_all_paths(
    nodes, edges, lines, linepaths, demands, demand_data, frequencies, create_env
):
    """Strategy 2:
    - multi commodity flow formulation for the passenger demand without any restrictions
    - multi-objective approach: first minimize cost then minimize travel times for passengers by allowing a cost
      degradation of 20%
    - add some cuts to the formulation to improve LP relaxation
    """
    logger.info("Starting line optimization using strategy 2.")
    with create_env() as env, gp.Model(env=env) as model:
        obj_cost = 0  # objective function for cost
        objtime = 0  # objective function for travel time

        # add variables for lines and frequencies
        x = {}
        numlines = 0
        for l in linepaths.index:
            numlines += 1
            for f in frequencies:
                x[l, f] = model.addVar(vtype=gp.GRB.BINARY, name=str(l) + str(f))
                obj_cost += x[l, f] * (
                    f * lines[l]["operating_cost"] + lines[l]["fix_cost"]
                )

        logger.info(
            f"Consider network with {len(nodes)} nodes, "
            f"{len(edges)} edges, and "
            f"{numlines} potential line routes"
        )

        # variable and objectve value for flow (flow from origin node s over edge (u,v))
        f = {}
        for s in nodes:
            for u, v in edges:
                f[s, u, v] = model.addVar(
                    lb=0,
                    vtype=gp.GRB.CONTINUOUS,
                    name="passflow" + str(s) + "," + str(u) + "," + str(v),
                )
                objtime += f[s, u, v] * edges[(u, v)]["time"]

        # each line with at most one frequency
        model.addConstrs(
            (gp.quicksum(x[i, j] for j in frequencies) <= 1 for i in linepaths.index),
            name="one-freq",
        )

        # add capacity constraints
        for u, v in edges:
            model.addConstr(
                (
                    gp.quicksum(
                        x[i, j] * lines[i]["capacity"] * j
                        for i in linepaths.index
                        for j in frequencies
                        if (u, v) in linepaths[i]
                    )
                )
                >= gp.quicksum(f[s, u, v] for s in nodes),
                name="capacity" + str(u) + "," + str(v),
            )

        # add constraints for multi commodity flow
        for s in nodes:
            # total demand from node s
            model.addConstr(
                gp.quicksum(f[s, s, u] for u in nodes if (s, u) in edges)
                == gp.quicksum(demands[s, v] for v in nodes if (s, v) in demands),
                name="demand" + str(s),
            )
            for v in nodes:
                if (s, v) in demands:
                    model.addConstr(
                        gp.quicksum(f[s, u, v] for u in nodes if (u, v) in edges)
                        == demands[s, v]
                        + gp.quicksum(f[s, v, w] for w in nodes if (v, w) in edges),
                        name="flowFrom-" + str(s) + ",to" + str(v),
                    )
                elif s != v:
                    model.addConstr(
                        gp.quicksum(f[s, u, v] for u in nodes if (u, v) in edges)
                        == gp.quicksum(f[s, v, w] for w in nodes if (v, w) in edges),
                        name="flowFrom-" + str(s) + ",to" + str(v),
                    )

        # Set up primary objective: minimize total cost
        model.setObjectiveN(obj_cost, index=0, priority=2, reltol=0.2, name="Cost")

        # Set up secondary objective: minimize passengers travel time
        model.setObjectiveN(objtime, index=1, priority=1, name="Traveltime")

        # add additional cuts to improve the LP relaxation: capacity and connectivity cuts around each station
        # where passengers want to start or end their trip
        for s in nodes:
            demand_from_s = demand_data.loc[demand_data["source"] == s, "demand"].sum()
            demand_to_s = demand_data.loc[demand_data["target"] == s, "demand"].sum()
            cap_from_s = 0
            lines_from_s = 0
            cap_to_s = 0
            lines_to_s = 0
            for l in linepaths.index:
                found_u = False
                found_v = False
                for u, v in linepaths[l]:
                    if s == u:
                        for f in frequencies:
                            cap_from_s += x[l, f] * lines[l]["capacity"] * f
                        if found_u == False:
                            for f in frequencies:
                                lines_from_s += x[l, f]
                            found_u = True
                    if s == v:
                        for f in frequencies:
                            cap_to_s += x[l, f] * lines[l]["capacity"] * f
                        if found_v == False:
                            for f in frequencies:
                                lines_to_s += x[l, f]
                            found_v = True
            if demand_from_s > 0:
                model.addConstr(cap_from_s >= demand_from_s, name="capFrom" + str(s))
                model.addConstr(lines_from_s >= 1, name="linesFrom" + str(s))
            if demand_to_s > 0:
                model.addConstr(cap_to_s >= demand_to_s, name="capTo" + str(s))
                model.addConstr(lines_to_s >= 1, name="linesTo" + str(s))

        model.optimize()

        # prepare return values
        obj_cost = -1
        lines_out = []
        if model.Status == gp.GRB.OPTIMAL:
            obj_cost = model.objVal
            # model.setParam(gp.GRB.Param.ObjNumber, 1)
            # objTime = model.objVal
            # save optimal line plan (after both objectives) in list of tuples
            # print non-zero solution values of line variables
            for i in linepaths.index:
                for j in frequencies:
                    if x[i, j].X > 0.5:
                        lines_out.append((i, j))
        return obj_cost, lines_out


def plot_lineplan(
    node_data: pd.DataFrame,
    edge_data: pd.DataFrame,
    linepath_data: pd.DataFrame,
    line_plan: list,
):
    """Visualize a line plan. The figure is opened in a browser
    A colormap with 20 different colors is used.
    If the line plan contains at most 10 lines, a colormap with 10 colors is used to have the colors more distinguishable.
    If the line plan contains more than 20 lines, the same colors are rerun.

    Parameters
    ----------
    node_data : DataFrame
        DataFrame with information on the nodes/stations. The frame must include "source".
    edge_data : DataFrame
        DataFrame with edges / connections between stations and associated attributes.
        The frame must include "source", "target", and "time"
    linepath_data : DataFrame
        DataFrame with information on the line routes/paths.
        It must include "linename", "edge_source", and "edge_target".
    line_plan: List
        A solution of the line optimization, i.e., a list with linenames and associated frequencies.

    """
    if mpl is None or nx is None:
        raise RuntimeError(
            "Plot not possible: networkx, matplotlib and matplotlib.pyplot are required for plotting the line plan"
        )

    if "posx" not in node_data.columns or "posy" not in node_data.columns:
        raise ValueError("Need posx and posy information in node_data!")

    if len(line_plan) > 20:
        logger.info(
            "Note that only 20 different colors are used. Line plan has more than 20 lines, hence, different lines will have the same color."
        )

    linepaths = (
        linepath_data.set_index("linename")
        .groupby(["linename"])
        .apply(lambda x: [(k, v) for k, v in zip(x["edge_source"], x["edge_target"])])
    )
    G = nx.from_pandas_edgelist(edge_data.reset_index(), create_using=nx.Graph())
    for number, row in node_data.set_index("number").iterrows():
        # print(number)
        G.add_node(number, pos=(row["posx"], row["posy"]))

    mpl.use("WebAgg")
    plt.figure(1, figsize=(14, 7))
    # plot network on the left
    plt.subplot(1, 2, 1)  # row 1, col 2 index 1
    pos = nx.get_node_attributes(G, "pos")
    nx.draw(G, pos, width=1, node_size=100, node_color="gray")

    if len(line_plan) <= 10:
        colormap = mpl.cm.tab10.colors
    else:
        colormap = mpl.cm.tab20.colors
    colornum = 0

    # plot line plan on the right
    plt.subplot(1, 2, 2)  # index 2
    plt.axis("off")
    # these parameters are used to "shift a line" when different lines use the same arc/edge
    # this is done in a very simple way; it might be the case that different lines
    # on the same edge are not easily distinguishable
    xmean = round(node_data.loc[:, "posx"].mean() / 50)
    ymean = round(node_data.loc[:, "posy"].mean() / 50)
    pathList = []
    for l, f in line_plan:
        lastU = -1
        for u, v in linepaths[l]:
            if lastU == v:
                break  # line path back
            lastU = u
            (x1, y1) = G.nodes[u]["pos"]
            (x2, y2) = G.nodes[v]["pos"]
            cnt = pathList.count((u, v)) + pathList.count((v, u))
            if abs(x1 - x2) < abs(y1 - y2):
                x1 += cnt * xmean
                x2 += cnt * xmean
            else:
                y1 += cnt * ymean
                y2 += cnt * ymean
            plt.plot((x1, x2), (y1, y2), linewidth=2, color=colormap[colornum])
            pathList.append((u, v))
        colornum += 1
        # start with first color if the number of colors is reached
        if colornum == len(colormap):
            colornum = 0

    # plot again all nodes
    for n in G.nodes:
        (x, y) = G.nodes[n]["pos"]
        plt.plot(
            x,
            y,
            marker="o",
            markersize=10,
            markeredgecolor="black",
            markerfacecolor="gray",
        )

    plt.show()
