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
    shortestPaths: bool = True,
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
        It must include "linename", "capacity", "fixCost", and "operatingCost".
    linepath_data : DataFrame
        DataFrame with information on the line routes/paths.
        It must include "linename", "edgeSource", and "edgeTarget".
    frequency: List
        List with possible frequencies: How often the line can be operated in the considered
        time horizon.

    Returns
    -------
    tuple
        - Cost of the optimal line concept (line frequency)
        - list of line-frequency tuples (optimal line concept)
    """
    missingData = False
    # check data
    if "number" not in node_data.columns:
        logger.info("column nr not present in node_data")
        missingData = True
    if "time" not in edge_data.columns:
        logger.info("column time not present in edge_data")
        missingData = True
    if "source" not in edge_data.columns:
        logger.info("column source not present in edge_data")
        missingData = True
    if "target" not in edge_data.columns:
        logger.info("column target not present in edge_data")
        missingData = True
    if "linename" not in line_data.columns:
        logger.info("column linename not present in line_data")
        missingData = True
    if "capacity" not in line_data.columns:
        logger.info("column capacity not present in line_data")
        missingData = True
    if "fixCost" not in line_data.columns:
        logger.info("column fixCost not present in line_data")
        missingData = True
    if "operatingCost" not in line_data.columns:
        logger.info("column operatingCost not present in line_data")
        missingData = True
    if "linename" not in linepath_data.columns:
        logger.info("column linename not present in linepath_data")
        missingData = True
    if "edgeSource" not in linepath_data.columns:
        logger.info("column edgeSource not present in linepath_data")
        missingData = True
    if "edgeTarget" not in linepath_data.columns:
        logger.info("column edgeTarget not present in linepath_data")
        missingData = True
    if "source" not in demand_data.columns:
        logger.info("column source not present in demand_data")
        missingData = True
    if "target" not in demand_data.columns:
        logger.info("column target not present in demand_data")
        missingData = True
    if "demand" not in demand_data.columns:
        logger.info("column demand not present in demand_data")
        missingData = True
    if node_data.isnull().values.any():
        logger.info("some value is nan or format is not correct for node_data")
        missingData = True
    if edge_data.isnull().values.any():
        logger.info("some value is nan or format is not correct for edge_data")
        missingData = True
    if line_data.isnull().values.any():
        logger.info("some value is nan or format is not correct for line_data")
        missingData = True
    if linepath_data.isnull().values.any():
        logger.info("some value is nan or format is not correct for linepath_data")
        missingData = True
    if demand_data.isnull().values.any():
        logger.info("some value is nan or format is not correct for demand_data")
        missingData = True
    if (demand_data < 0).values.any():
        logger.info(
            "All demand should be non-negative! Some value is negative in demand_data"
        )
        missingData = True
    if missingData:
        return -1, []

    # prepare data
    nodes = node_data.set_index("number").to_dict("index")
    edges = edge_data.set_index(["source", "target"]).to_dict("index")
    lines = line_data.set_index("linename").to_dict("index")
    linepaths = (
        linepath_data.set_index("linename")
        .groupby(["linename"])
        .apply(lambda x: [(k, v) for k, v in zip(x["edgeSource"], x["edgeTarget"])])
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

    if shortestPaths:
        if nx is None:
            logger.info(
                "Networkx is needed for strategy 1 but not available. Using strategy 2 instead."
            )
            return allowAllPaths(
                nodes,
                edges,
                lines,
                linepaths,
                demands,
                demand_data,
                frequencies,
                create_env,
            )
        return allShortestPaths(
            nodes, edges, edge_data, lines, linepaths, demands, frequencies, create_env
        )
    else:
        return allowAllPaths(
            nodes,
            edges,
            lines,
            linepaths,
            demands,
            demand_data,
            frequencies,
            create_env,
        )


# Strategy 1:
#  - compute all shortest path for each OD pair using networkx all_shortest_paths algorithm
#  - create passenger variable for each path
def allShortestPaths(
    nodes, edges, edge_data, lines, linepaths, demands, frequencies, create_env
):
    logger.info("Starting line optimization using strategy 2.")
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
                    obj=f * lines[l]["operatingCost"] + lines[l]["fixCost"],
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
        f = {}  # stuv from s to t using edge uv
        for s, t in demands:
            paths = list(nx.all_shortest_paths(G, source=s, target=t, weight="time"))
            demandExpr = 0
            for p in paths:
                y = model.addVar(
                    ub=demands[s, t],
                    vtype=gp.GRB.CONTINUOUS,
                    name="passpath" + str(s) + "," + str(t) + str(p),
                )
                demandExpr += y
                for node in range(len(p) - 1):
                    model.chgCoeff(cap[p[node], p[node + 1]], y, -1.0)
            model.addConstr(demandExpr == demands[s, t])

        # each line with at most one frequency
        model.addConstrs(
            (gp.quicksum(x[i, j] for j in frequencies) <= 1 for i in linepaths.index),
            name="one-freq",
        )

        model.optimize()

        # prepare return values
        objCost = -1
        linesOut = []
        if model.Status == gp.GRB.OPTIMAL:
            objCost = model.objVal
            for i in linepaths.index:
                for j in frequencies:
                    if x[i, j].X > 0.5:
                        linesOut.append((i, j))

        return objCost, linesOut


# Strategy 2:
#  - multi commodity flow formulation for the passenger demand without any restrictions
#  - multi-objective approach: first minimize cost then minimize travel times for passengers by allowing a cost
#    degradation of 20%
#  - add some cuts to the formulation to improve LP relaxation
def allowAllPaths(
    nodes, edges, lines, linepaths, demands, demand_data, frequencies, create_env
):
    logger.info("Starting line optimization using strategy 1.")
    with create_env() as env, gp.Model(env=env) as model:
        objcost = 0  # objective function for cost
        objtime = 0  # objective function for travel time

        # add variables for lines and frequencies
        x = {}
        numlines = 0
        for l in linepaths.index:
            numlines += 1
            for f in frequencies:
                x[l, f] = model.addVar(vtype=gp.GRB.BINARY, name=str(l) + str(f))
                objcost += x[l, f] * (
                    f * lines[l]["operatingCost"] + lines[l]["fixCost"]
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
        model.setObjectiveN(objcost, index=0, priority=2, reltol=0.2, name="Cost")

        # Set up secondary objective: minimize passengers travel time
        model.setObjectiveN(objtime, index=1, priority=1, name="Traveltime")

        # add additional cuts to improve the LP relaxation: capacity and connectivity cuts around each station
        # where passengers want to start or end their trip
        for s in nodes:
            demandFromS = demand_data.loc[demand_data["source"] == s, "demand"].sum()
            demandToS = demand_data.loc[demand_data["target"] == s, "demand"].sum()
            capFS = 0
            linesFS = 0
            capTS = 0
            linesTS = 0
            for l in linepaths.index:
                foundU = False
                foundV = False
                for u, v in linepaths[l]:
                    if s == u:
                        for f in frequencies:
                            capFS += x[l, f] * lines[l]["capacity"] * f
                        if foundU == False:
                            for f in frequencies:
                                linesFS += x[l, f]
                            foundU = True
                    if s == v:
                        for f in frequencies:
                            capTS += x[l, f] * lines[l]["capacity"] * f
                        if foundV == False:
                            for f in frequencies:
                                linesTS += x[l, f]
                            foundV = True
            if demandFromS > 0:
                model.addConstr(capFS >= demandFromS, name="capFrom" + str(s))
                model.addConstr(linesFS >= 1, name="linesFrom" + str(s))
            if demandToS > 0:
                model.addConstr(capTS >= demandToS, name="capTo" + str(s))
                model.addConstr(linesTS >= 1, name="linesTo" + str(s))

        model.optimize()

        # prepare return values
        objCost = -1
        linesOut = []
        if model.Status == gp.GRB.OPTIMAL:
            objCost = model.objVal
            # model.setParam(gp.GRB.Param.ObjNumber, 1)
            # objTime = model.objVal
            # save optimal line plan (after both objectives) in list of tuples
            # print non-zero solution values of line variables
            for i in linepaths.index:
                for j in frequencies:
                    if x[i, j].X > 0.5:
                        linesOut.append((i, j))
        return objCost, linesOut
