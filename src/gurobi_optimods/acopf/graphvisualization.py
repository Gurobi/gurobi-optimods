# https://networkx.org/documentation/stable/_downloads/networkx_reference.pdf

import plotly.graph_objects as go
import plotutils as pu
import psutil 
import numpy as np
import networkx as nx
from collections import defaultdict
import plotly.graph_objects as go
from typing import Any, List, Dict, Tuple, Union, Callable

Vertex = Any
Edge = Tuple[Vertex, Vertex]
Num = Union[int, float]

class GraphVisualization:
    """Description"""#FIXME add more details

    def __init__(
        self,
        G: nx.Graph,
        pos: Dict[Vertex, Union[Tuple[Num, Num], Tuple[Num, Num, Num]]],
        node_text: Union[Dict[Vertex, str], Callable] = None,
        node_text_position: Union[Dict[Vertex, str], Callable, str] = None,
        node_text_font_color: Union[Dict[Vertex, str], Callable, str] = None,
        node_text_font_family: Union[Dict[Vertex, str], Callable, str] = None,
        node_text_font_size: Union[Dict[Vertex, Num], Callable, str] = None,
        node_size: Union[Dict[Vertex, Num], Callable, Num] = None,
        node_color: Union[Dict[Vertex, Union[str, Num]], Callable, Union[str, Num]] = None,
        node_border_width: Union[Dict[Vertex, Num], Callable, Num] = None,
        node_border_color: Union[Dict[Vertex, str], Callable, str] = None,
        node_opacity: Num = None,
        edge_width: Union[Dict[Edge, Num], Callable, Num] = None,
        edge_color: Union[Dict[Edge, str], Callable, str] = None,
        edge_opacity: Num = None,
    ):
        # check dimensions
        if all(len(pos.get(v, [])) == 2 for v in G):
            self.is_3d = False
        elif all(len(pos.get(v, [])) == 3 for v in G):
            self.is_3d = True
        else:
            raise ValueError

        # default settings
        self.default_settings = dict(
            node_text=str,  # show node label
            node_text_position="middle center",
            node_text_font_color='#000000',
            node_text_font_family='Arial',
            node_text_font_size=14,
            node_size=10 if self.is_3d else 18,
            node_color='#fcfcfc',
            node_border_width=2,
            node_border_color='#333333',
            node_opacity=0.8,
            edge_width=4 if self.is_3d else 2,
            edge_color='#000000', #'#808080',
            edge_opacity=0.8,
        )

        # save settings
        self.G = G
        self.pos = pos
        self.node_text = node_text
        self.node_text_position = node_text_position
        self.node_text_font_color = node_text_font_color
        self.node_text_font_family = node_text_font_family
        self.node_text_font_size = node_text_font_size
        self.node_size = node_size
        self.node_color = node_color
        self.node_border_width = node_border_width
        self.node_border_color = node_border_color
        self.node_opacity = node_opacity
        self.edge_width = edge_width
        self.edge_color = edge_color
        self.edge_opacity = edge_opacity

    def _get_edge_traces(self) -> List[Union[go.Scatter, go.Scatter3d]]:
        # group all edges by (color, width)
        groups = defaultdict(list)

        for edge in self.G.edges():
            color = self._get_setting('edge_color', edge)
            width = self._get_setting('edge_width', edge)
            groups[(color, width)] += [edge]

        # process each group
        traces = []
        for (color, width), edges in groups.items():
            x, y, z = [], [], []
            for v, u in edges:
                x += [self.pos[v][0], self.pos[u][0], None]
                y += [self.pos[v][1], self.pos[u][1], None]
                if self.is_3d:
                    z += [self.pos[v][2], self.pos[u][2], None]

            params = dict(
                x=x,
                y=y,
                mode='lines',
                hoverinfo='none',
                line=dict(color=color, width=width),
                opacity=self._get_setting('edge_opacity'),
            )

            traces += [go.Scatter3d(z=z, **params) if self.is_3d else go.Scatter(**params)]

        return traces

    def _get_node_trace(self, showlabel, colorscale, showscale, colorbar_title, reversescale) -> Union[go.Scatter, go.Scatter3d]:
        x, y, z = [], [], []
        for v in self.G.nodes():
            x += [self.pos[v][0]]
            y += [self.pos[v][1]]
            if self.is_3d:
                z += [self.pos[v][2]]

        params = dict(
            x=x,
            y=y,
            mode='markers' + ('+text' if showlabel else ''),
            hoverinfo='text',
            marker=dict(
                showscale=showscale,
                colorscale=colorscale,
                reversescale=reversescale,
                color=self._get_setting('node_color'),
                size=self._get_setting('node_size'),
                line_width=self._get_setting('node_border_width'),
                line_color=self._get_setting('node_border_color'),
                colorbar=dict(
                    thickness=15,
                    title=colorbar_title,
                    xanchor='left',
                    titleside='right'
                ),
            ),
            text=self._get_setting('node_text'),
            textfont=dict(
                color=self._get_setting('node_text_font_color'),
                family=self._get_setting('node_text_font_family'),
                size=self._get_setting('node_text_font_size')
            ),
            textposition=self._get_setting('node_text_position'),
            opacity=self._get_setting('node_opacity'),
        )

        trace = go.Scatter3d(z=z, **params) if self.is_3d else go.Scatter(**params)
        return trace

    def _get_setting(self, setting_name, edge=None):
        default_setting = self.default_settings.get(setting_name)
        def_func = default_setting if callable(default_setting) else lambda x: default_setting
        setting = self.__dict__.get(setting_name)

        if edge is None:  # vertex-specific
            if setting is None:  # default is used
                if callable(default_setting):  # default is a function
                    return [def_func(v) for v in self.G.nodes()]
                else:  # default is a constant
                    return default_setting
            elif callable(setting):  # setting is a function
                return [setting(v) for v in self.G.nodes()]
            elif isinstance(setting, dict):  # setting is a dict
                return [setting.get(v, def_func(v)) for v in self.G.nodes()]
            else:  # setting is a constant
                return setting
        else:  # edge-specific
            if setting is None:  # default is used
                return def_func(edge)
            elif callable(setting):  # setting is a function
                return setting(edge)
            elif isinstance(setting, dict):  # setting is a dict
                return setting.get(edge, def_func(edge))
            else:  # setting is a constant
                return setting

    def create_figure(
        self,
        showlabel=True,
        colorscale='YlGnBu',
        showscale=False,
        colorbar_title='',
        reversescale=False,
        **params
    ) -> go.Figure:
        axis_settings = dict(
            autorange=True,
            showgrid=False,
            zeroline=False,
            showline=False,
            visible=False,
            ticks='',
            showticklabels=False,
        )
        scene = dict(
            xaxis=axis_settings,
            yaxis=axis_settings,
            zaxis=axis_settings,
        )

        layout_params = dict(
            paper_bgcolor='rgba(255,255,255,255)',  # white
            plot_bgcolor='rgba(0,0,0,0)',  # transparent
            autosize=False,
            height=400,
            width=450 if showscale else 375,
            title='',
            titlefont_size=16,
            showlegend=False,
            hovermode='closest',
            margin=dict(b=5, l=0, r=0, t=20),
            annotations=[],
            xaxis=axis_settings,
            yaxis=axis_settings,
            scene=scene,
        )

        # override with the given parameters
        layout_params.update(params)

        # create figure
        fig = go.Figure(layout=go.Layout(**layout_params))
        fig.add_traces(self._get_edge_traces())
        fig.add_trace(self._get_node_trace(showlabel, colorscale, showscale, colorbar_title, reversescale))
        return fig
