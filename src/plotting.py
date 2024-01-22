import matplotlib.pyplot as plt
import numpy as np
from src.data_collections import Geometry, PanelizedGeometry
import typing as tp


def plot_panelized_geometry(geometry: Geometry, panelized_geometry: PanelizedGeometry):
    fig = plt.figure(1)  # Create figure
    plt.cla()
    plt.fill(geometry.x, geometry.y, 'k')  # Plot polygon (circle or airfoil)
    X = panelized_geometry.control_points_x_cor + panelized_geometry.panel_length * np.cos(
        panelized_geometry.panel_normal_angle)
    Y = panelized_geometry.control_points_y_cor + panelized_geometry.panel_length * np.sin(
        panelized_geometry.panel_normal_angle)
    number_of_panels = len(X)  # Number of panels
    for i in range(number_of_panels):
        if (i == 0):  # For first panel
            plt.plot([panelized_geometry.control_points_x_cor[i], X[i]],
                     [panelized_geometry.control_points_y_cor[i], Y[i]],
                     'b-', label='First Panel')  # Plot the first panel normal vector
        elif (i == 1):  # For second panel
            plt.plot([panelized_geometry.control_points_x_cor[i], X[i]],
                     [panelized_geometry.control_points_y_cor[i], Y[i]],
                     'g-', label='Second Panel')  # Plot the second panel normal vector
        else:  # For every other panel
            plt.plot([panelized_geometry.control_points_x_cor[i], X[i]],
                     [panelized_geometry.control_points_y_cor[i], Y[i]],
                     'r-')

    plt.xlabel('X-Axis')  # Set X-label
    plt.ylabel('Y-Axis')  # Set Y-label
    plt.title('Panel Geometry')  # Set title
    plt.axis('equal')  # Set axes equal
    plt.legend()  # Plot legend
    return fig


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
