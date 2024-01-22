import numpy as np
import matplotlib.pyplot as plt
import typing as tp
import multiprocessing as mp
from itertools import repeat
from . import elementary_flows
from .plotting import plot_flow_from_stream_function, plot_flow_from_velocities

plotting_kwargs = {
    'X_NEG_LIMIT': -5,
    'X_POS_LIMIT': 5,
    'Y_NEG_LIMIT': -5,
    'Y_POS_LIMIT': 5,
    'STREAMLINE_DENSITY': 3,
    'STREAMLINE_COLOR': 'b',
    'CONTOR_LEVELS': 50,
    'CONTOR_COLOR': 'k',
    'FIGURE_SIZE': (12, 12),
    'DPI': 100,
    "CONTOUR_LABELS": True
}


def velocity_internal(flow, x, y):
    return flow.velocity(x, y)


def stream_function_internal(flow, x, y):
    return flow.stream_function(x, y)


class FlowField:
    '''
    A class that represents a flow field.

    Attributes
    ----------

    flows : list
    A list of elementary flows that make up the flow field.

    Methods
    -------

    stream_function(x, y)
    Returns the stream function of the flow field at the point (x, y).

    velocity(x, y)
    Returns the velocity of the flow field at the point (x, y).

    plot()
    Plots the stream function of the flow field.

    plot_velocity()
    Plots the velocity of the flow field as streamlines.

    '''

    def __init__(self, flows: tp.Optional[tp.List[elementary_flows.ElementaryFlow]] = None, **kwargs):
        if flows is None:
            flows = []
        self.flows: tp.Optional[tp.List[elementary_flows.ElementaryFlow]] = flows
        self.plotting_kwargs: dict = plotting_kwargs

        self.plotting_kwargs.update(**kwargs)

    def stream_function(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        pool = mp.Pool(mp.cpu_count())
        return sum(pool.starmap(stream_function_internal, zip(self.flows, repeat(x), repeat(y))))
        # return sum([flow.stream_function(x, y) for flow in self.flows])

    def plot_flow_from_stream_function(self, x: np.ndarray, y: np.ndarray) -> plt.Figure:
        X, Y = np.meshgrid(x, y)
        return plot_flow_from_stream_function(self.stream_function, X, Y, **self.plotting_kwargs)

    def velocity(self, x: np.ndarray, y: np.ndarray) -> tp.Tuple[np.ndarray, np.ndarray]:
        pool = mp.Pool(mp.cpu_count())

        flow_velocities = pool.starmap(velocity_internal, zip(self.flows, repeat(x), repeat(y)))
        # flow_velocities = [flow.velocity(x, y) for flow in self.flows]
        U = sum([flow_vel[0] for flow_vel in flow_velocities])
        V = sum([flow_vel[1] for flow_vel in flow_velocities])
        return U, V

    def plot_velocity(self, x: np.ndarray, y: np.ndarray) -> plt.Figure:
        X, Y = np.meshgrid(x, y)
        U, V = self.velocity(X, Y)
        return plot_flow_from_velocities(X, Y, U, V, **self.plotting_kwargs)
