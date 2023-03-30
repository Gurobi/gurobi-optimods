import plotly.graph_objects as go
import plotutils as pu
import psutil
import logging
import numpy as np
from collections import defaultdict
import plotly.graph_objects as go
from typing import Any, List, Dict, Tuple, Union, Callable

from .grbgraph import *
from .myutils import break_exit

Vertex = Any
Edge = Tuple[Vertex, Vertex]
Num = Union[int, float]


class plotlyhandler:
    """Description"""

    #
    #
    # Creates a plotly graph to be rendered
    #
    # Input is a graph given by coordinates
    #
    def __init__(
        self,
        gG: grbGraph,
        pos,
        annotation_list,
        vertex_size=None,
        vertex_color=None,
        edge_width=None,
        edge_color=None,
        edge_map=None,
        vertex_text=None,
        vertex_text_position=None,
        vertex_text_font_color=None,
        vertex_text_font_family=None,
        vertex_text_font_size=None,
        vertex_border_width=None,
        vertex_border_color=None,
        vertex_opacity=None,
        edge_opacity=None,
    ):

        if vertex_opacity == None:
            vertex_opacity = 0.8
        if edge_opacity == None:
            edge_opacity = 0.8
        if vertex_text_position == None:
            vertex_text_position = "middle center"
        if vertex_text_font_color == None:
            vertex_text_font_color = ("#000000",)
        if vertex_text_font_family == None:
            vertex_text_font_family = "Arial"
        if vertex_text_font_size == None:
            vertex_text_font_size = 14

        self.gG = gG
        self.pos = pos
        self.annotation_list = annotation_list
        self.vertex_text_position = vertex_text_position
        self.vertex_text_font_color = vertex_text_font_color
        self.vertex_text_font_family = vertex_text_font_family
        self.vertex_text_font_size = vertex_text_font_size
        self.vertex_text = vertex_text
        self.vertex_size = vertex_size
        self.vertex_color = vertex_color
        self.vertex_border_width = vertex_border_width
        self.vertex_border_color = vertex_border_color
        self.vertex_opacity = vertex_opacity
        self.edge_width = edge_width
        self.edge_color = edge_color
        self.input_edge_width = edge_width
        self.input_edge_color = edge_color
        self.edge_opacity = edge_opacity
        self.edge_map = edge_map

    def generate_edge_traces(self):  # -> List[Union[go.Scatter, go.Scatter3d]]:
        # group all edges by (color, width)
        groups = defaultdict(list)

        # print(self.input_edge_color)

        count = 1
        localdeg = {}
        # for edge in self.G.edges():
        for edgecnt in range(self.gG.m):
            edge = self.gG.edges[edgecnt]
            small = min(edge[0] + 1, edge[1] + 1)
            large = max(edge[0] + 1, edge[1] + 1)
            localdeg[small, large] = 0

        loud = False

        # for edge in self.G.edges():
        for edgecnt in range(self.gG.m):
            edge = self.gG.edges[edgecnt]
            small = min(edge[0] + 1, edge[1] + 1)
            large = max(edge[0] + 1, edge[1] + 1)
            thedeg = localdeg[small, large]
            position = self.edge_map[(small, large, thedeg)]

            width = self.input_edge_width[position]
            color = self.input_edge_color[position]

            if loud:

                logging.info(
                    "edge %d is (+1) (%d, %d) position %d width %d localdeg %d\n"
                    % (
                        count,
                        edge[0] + 1,
                        edge[1] + 1,
                        position,
                        self.input_edge_width[position],
                        localdeg[small, large],
                    )
                )
                # print('edge',count, 'is (+1) (',edge[0]+1,',',edge[1]+1, ') position', position, 'w', self.input_edge_width[position], 'localdeg', localdeg[small,large])

            localdeg[small, large] += 1

            groups[(color, width)] += [edge]

            # print('edge',count, 'is (+1)',edge[0]+1,',',edge[1]+1, 'of width', width)
            count += 1

        # process each group
        traces = []
        for (color, width), edges in groups.items():
            x, y, z = [], [], []
            for v, u in edges:
                x += [self.pos[v][0], self.pos[u][0], None]
                y += [self.pos[v][1], self.pos[u][1], None]

            params = dict(
                x=x,
                y=y,
                mode="lines",
                hoverinfo="none",
                line=dict(color=color, width=width),
                opacity=self.edge_opacity,
            )

            traces += [go.Scatter(**params)]

        return traces

    def generate_vertex_trace(
        self, showlabel, colorscale, showscale, colorbar_title, reversescale
    ) -> Union[go.Scatter, go.Scatter3d]:
        x, y, z = [], [], []
        for v in self.gG.vertices.values():
            x += [self.pos[v][0]]
            y += [self.pos[v][1]]

        params = dict(
            x=x,
            y=y,
            mode="markers" + ("+text" if showlabel else ""),
            hoverinfo="text",
            marker=dict(
                showscale=showscale,
                colorscale=colorscale,
                reversescale=reversescale,
                color=self.vertex_color_list(),
                size=self.vertex_size_list(),
                line_width=self.vertex_border_width_list(),
                line_color=self.vertex_border_color,
                colorbar=dict(
                    thickness=15,
                    title=colorbar_title,
                    xanchor="left",
                    titleside="right",
                ),
            ),
            text=self.vertex_text_list(),
            textfont=dict(
                color=self.vertex_text_font_color,
                family=self.vertex_text_font_family,
                size=self.vertex_text_font_size,
            ),
            textposition=self.vertex_text_position,
            opacity=self.vertex_opacity,
        )

        trace = go.Scatter(params)
        return trace

    def vertex_color_list(self):
        return [self.vertex_color[v] for v in self.gG.vertices.values()]

    def vertex_size_list(self):
        return [self.vertex_size[v] for v in self.gG.vertices.values()]

    def vertex_border_width_list(self):
        return [self.vertex_border_width[v] for v in self.gG.vertices.values()]

    def vertex_text_list(self):
        return [self.vertex_text[v] for v in self.gG.vertices.values()]  #

    def create_figure(
        self,
        showlabel=True,
        colorscale="YlGnBu",
        showscale=False,
        colorbar_title="",
        reversescale=False,
        **params
    ) -> go.Figure:
        axis_settings = dict(
            autorange=True,
            showgrid=False,
            zeroline=False,
            showline=False,
            visible=False,
            ticks="",
            showticklabels=False,
        )
        scene = dict(
            xaxis=axis_settings,
            yaxis=axis_settings,
            zaxis=axis_settings,
        )

        layout_params = dict(
            paper_bgcolor="rgba(255,255,255,255)",
            plot_bgcolor="rgba(0,0,0,0)",
            autosize=False,
            height=400,
            width=450,
            title="",
            titlefont_size=16,
            showlegend=False,
            hovermode="closest",
            margin=dict(b=5, l=0, r=0, t=20),
            annotations=[],
            xaxis=axis_settings,
            yaxis=axis_settings,
            scene=scene,
        )

        layout_params.update(params)

        if False:
            print([self.vertex_size[v] for v in self.gG.vertices.values()])
            print([self.vertex_color[v] for v in self.gG.vertices.values()])
            break_exit("prefig")

        # create figure and add traces
        fig = go.Figure(layout=go.Layout(**layout_params))
        fig.add_traces(self.generate_edge_traces())
        fig.add_trace(
            self.generate_vertex_trace(
                showlabel, colorscale, showscale, colorbar_title, reversescale
            )
        )
        if len(self.annotation_list) > 0:
            for k in range(len(self.annotation_list)):
                fig.add_annotation(
                    x=200,
                    y=2000 - 70 * k,
                    text=self.annotation_list[k],  # "BALONEY",
                    showarrow=False,
                    font=dict(size=20, color="red"),
                    yshift=10,
                )

        return fig
