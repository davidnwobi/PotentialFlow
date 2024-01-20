import numpy as np
import matplotlib.pyplot as plt
import elementary_flows
import typing as tp

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


def plot_flow_from_stream_function(psi: tp.Callable[[np.ndarray, np.ndarray], np.ndarray], X: np.ndarray,
                                   Y: np.ndarray, **kwargs) -> plt.Figure:
    fig = plt.figure(figsize=kwargs.get("FIGURE_SIZE", (12, 12)), dpi=kwargs.get("DPI", 100))
    CS = plt.contour(X, Y, psi(X, Y), kwargs.get("CONTOR_LEVELS", 50))
    if kwargs.get("CONTOUR_LABELS"):
        plt.clabel(CS, inline=1, fontsize=10)
    plt.xlim(kwargs.get("X_NEG_LIMIT"), kwargs.get("X_POS_LIMIT"))
    plt.ylim(kwargs.get("Y_NEG_LIMIT"), kwargs.get("Y_POS_LIMIT"))
    return fig


def plot_flow_from_velocities(X: np.ndarray, Y: np.ndarray, U: np.ndarray, V: np.ndarray, **kwargs) -> plt.Figure:
    fig = plt.figure(figsize=kwargs.get("FIGURE_SIZE", (12, 12)), dpi=kwargs.get("DPI", 100))
    plt.streamplot(X, Y, U, V, density=kwargs.get("STREAMLINE_DENSITY", 3), color=kwargs.get("STREAMLINE_COLOR", 'b'))
    plt.xlim(kwargs.get("X_NEG_LIMIT"), kwargs.get("X_POS_LIMIT"))
    plt.ylim(kwargs.get("Y_NEG_LIMIT"), kwargs.get("Y_POS_LIMIT"))
    return fig


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
        self.plotting_kwargs = plotting_kwargs

        self.plotting_kwargs.update(**kwargs)

    def stream_function(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        return sum([flow.stream_function(x, y) for flow in self.flows])

    def plot_flow_from_stream_function(self, x: np.ndarray, y: np.ndarray) -> plt.Figure:
        X, Y = np.meshgrid(x, y)
        return plot_flow_from_stream_function(self.stream_function, X, Y, **self.plotting_kwargs)

    def velocity(self, x: np.ndarray, y: np.ndarray) -> tp.Tuple[np.ndarray, np.ndarray]:
        flow_velocities = [flow.velocity(x, y) for flow in self.flows]
        U = sum([flow_vel[0] for flow_vel in flow_velocities])
        V = sum([flow_vel[1] for flow_vel in flow_velocities])
        return U, V

    def plot_velocity(self, x: np.ndarray, y: np.ndarray) -> plt.Figure:
        X, Y = np.meshgrid(x, y)
        U, V = self.velocity(X, Y)
        return plot_flow_from_velocities(X, Y, U, V, **self.plotting_kwargs)

