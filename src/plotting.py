import matplotlib.pyplot as plt
import numpy as np
from src.data_collections import Geometry, PanelizedGeometry
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
            plt.plot([panelized_geometry.control_points_x_cor[i], X[i]], [panelized_geometry.control_points_y_cor[i], Y[i]],
                     'b-', label='First Panel')  # Plot the first panel normal vector
        elif (i == 1):  # For second panel
            plt.plot([panelized_geometry.control_points_x_cor[i], X[i]], [panelized_geometry.control_points_y_cor[i], Y[i]],
                     'g-', label='Second Panel')  # Plot the second panel normal vector
        else:  # For every other panel
            plt.plot([panelized_geometry.control_points_x_cor[i], X[i]], [panelized_geometry.control_points_y_cor[i], Y[i]],
                     'r-')

    plt.xlabel('X-Axis')  # Set X-label
    plt.ylabel('Y-Axis')  # Set Y-label
    plt.title('Panel Geometry')  # Set title
    plt.axis('equal')  # Set axes equal
    plt.legend()  # Plot legend
    return fig