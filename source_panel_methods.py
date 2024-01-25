import matplotlib.pyplot as plt
from src.panel_generator import PanelGenerator
from src.geometric_integrals import compute_panel_geometric_integrals_source
import src.data_collections as dc
from src.flow_field import *
from src.source_panel_methods_funcs import *
import numpy as np
import pandas as pd
from multiprocessing import freeze_support
import time

if __name__ == '__main__':
    time_start = time.time()
    freeze_support()
    numB = 500  # Number of boundary points
    num_grid = 1000
    X_NEG_LIMIT = -2
    X_POS_LIMIT = 2
    Y_NEG_LIMIT = -2
    Y_POS_LIMIT = 2
    tO = 22.5  # Angle offset [deg]
    load = 'Circle'  # Load circle or airfoil
    radius = 1
    theta = np.linspace(0, 360, numB)  # Angles to compute boundary points [deg]
    theta = theta + tO  # Add angle offset [deg]
    theta = theta * (np.pi / 180)  # Convert angle to radians [rad]

    # Boundary points
    XB = np.cos(theta) * radius  # X value of boundary points
    YB = np.sin(theta) * radius  # Y value of boundary points

    AoA = 20
    V = 1
    geometry = dc.Geometry(XB, YB, AoA)
    panelized_geometry = PanelGenerator.compute_geometric_quantities(geometry)

    panel_time_start = time.time()
    I, J = compute_panel_geometric_integrals_source(panelized_geometry)
    panel_time_end = time.time()
    print('Time to compute panels: ', panel_time_end - panel_time_start)

    gamma_time_start = time.time()
    gamma = compute_source_strengths(panelized_geometry, V, I)
    gamma_time_end = time.time()
    print('Time to compute gamma: ', gamma_time_end - gamma_time_start)
    sumLambda = sum(gamma * panelized_geometry.S)  # Check sum of source panel strengths
    print("Sum of L: ", sumLambda)

    V_normal, V_tangential = compute_panel_velocities_source(panelized_geometry, gamma, V, I, J)
    x, y = np.linspace(X_NEG_LIMIT, X_POS_LIMIT, num_grid), np.linspace(Y_NEG_LIMIT, Y_POS_LIMIT, num_grid)  # Create grid

    velocity_time_start = time.time()
    u, v = compute_grid_velocity_source(panelized_geometry, x, y, gamma, V, AoA)
    velocity_time_end = time.time()
    print('Time to compute grid velocity: ', velocity_time_end - velocity_time_start)

    X = panelized_geometry.xC - panelized_geometry.S / 2 * np.cos(panelized_geometry.phi)
    Y = panelized_geometry.yC - panelized_geometry.S / 2 * np.sin(panelized_geometry.phi)
    X = np.append(X, X[0])
    Y = np.append(Y, Y[0])
    panel_vector_X = panelized_geometry.xC + panelized_geometry.S / 2 * np.cos(panelized_geometry.delta)
    panel_vector_Y = panelized_geometry.yC + panelized_geometry.S / 2 * np.sin(panelized_geometry.delta)

    local_v = np.sqrt(u ** 2 + v ** 2)
    cp = 1 - (local_v / V) ** 2

    # Cp graph
    panel_velocities = V_tangential ** 2

    theoretical_Cp = lambda angle: 1 - 4 * (np.sin(angle)) ** 2
    theta = np.linspace(0, 2 * np.pi, 100)

    Cp = 1 - panel_velocities / V ** 2
    panelized_geometry.beta = np.where(panelized_geometry.beta > 2 * np.pi, panelized_geometry.beta - 2 * np.pi,
                                       panelized_geometry.beta)
    panelized_geometry.beta = np.where(panelized_geometry.beta < 0, panelized_geometry.beta + 2 * np.pi,
                                       panelized_geometry.beta)

    CN = -Cp * np.sin(panelized_geometry.beta) * panelized_geometry.S  # Normal coefficient
    CA = -Cp * np.cos(panelized_geometry.beta) * panelized_geometry.S  # Axial coefficient

    CL = np.sum(CN * np.cos(AoA * np.pi / 180)) - np.sum(CA * np.sin(AoA * np.pi / 180))  # Lift coefficient
    CD = np.sum(CN * np.sin(AoA * np.pi / 180)) + np.sum(CA * np.cos(AoA * np.pi / 180))  # Drag coefficient

    fig, axs = plt.subplots(2, 3, figsize=(30, 20), dpi=200)
    fig.delaxes(axs[1, 2])
    # Panel Geometry
    axs[0, 0].set_title('Panel Geometry')
    for i in range(len(panelized_geometry.S)):
        axs[0, 0].plot([X[i], X[i + 1]], [Y[i], Y[i + 1]], 'k')
        axs[0, 0].plot([X[i]], [Y[i]], 'ro')
        axs[0, 0].plot([panelized_geometry.xC[i]], [panelized_geometry.yC[i]], 'bo')
        axs[0, 0].text(panelized_geometry.xC[i], panelized_geometry.yC[i], str(i))
        if i == 0:
            axs[0, 0].plot([panelized_geometry.xC[i], panel_vector_X[i]], [panelized_geometry.yC[i], panel_vector_Y[i]],
                           'k', label='First Panel')
        if i == 1:
            axs[0, 0].plot([panelized_geometry.xC[i], panel_vector_X[i]], [panelized_geometry.yC[i], panel_vector_Y[i]],
                           'k', label='Second Panel')
        else:
            axs[0, 0].plot([panelized_geometry.xC[i], panel_vector_X[i]], [panelized_geometry.yC[i], panel_vector_Y[i]], 'k')
        axs[0, 0].plot([panel_vector_X[i]], [panel_vector_Y[i]], 'go')
    axs[0, 0].set_xlabel('X')
    axs[0, 0].set_ylabel('Y')
    axs[0, 0].legend()
    axs[0, 0].set_aspect('equal', adjustable='box')

    # Streamlines
    axs[0, 1].set_title('Streamlines')
    axs[0, 1].streamplot(x, y, u, v, density=2, color='b')
    axs[0, 1].set_xlim(X_NEG_LIMIT, X_POS_LIMIT)
    axs[0, 1].set_ylim(Y_NEG_LIMIT, Y_POS_LIMIT)
    axs[0, 1].fill(geometry.x, geometry.y, 'k')  # Plot polygon (circle or airfoil)
    axs[0, 1].set_xlabel('X')
    axs[0, 1].set_ylabel('Y')
    axs[0, 1].set_aspect('equal', adjustable='box')

    # Pressure Contours
    axs[0, 2].set_title('Pressure Contours')
    cp_plot = axs[0, 2].contourf(x, y, cp, 100, cmap='jet')
    plt.colorbar(cp_plot, ax=axs[0, 2], label="Pressure Coefficient")
    axs[0, 2].fill(geometry.x, geometry.y, 'k')
    axs[0, 2].set_xlabel('X')
    axs[0, 2].set_ylabel('Y')
    axs[0, 2].set_aspect('equal', adjustable='box')

    # Cp Distribution
    axs[1, 0].set_title('Cp Distribution')
    axs[1, 0].plot(theta, theoretical_Cp(theta), 'k', label='Theoretical')
    axs[1, 0].plot(panelized_geometry.beta, Cp, 'ro', label='Panel Method')
    axs[1, 0].set_xlabel('Theta [rad]')
    axs[1, 0].set_ylabel('Cp')
    axs[1, 0].legend()

    # Results Table
    axs[1, 1].axis('off')
    table_data = pd.DataFrame({
        "CL": [round(CL, 6)],
        "CD": [round(CD, 6)],
        "Sum of Source Strengths": [round(sumLambda, 6)]
    })
    table_data = table_data.round(6)
    table_data = table_data.T
    row_labels = table_data.index.values
    # Customize the table appearance
    table = axs[1, 1].table(cellText=table_data.values,
                            rowLabels=row_labels,
                            loc='center',
                            cellLoc='center',
                            bbox=[0.3, 0.6, 0.1, 0.2])

    # Adjust font size and style
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.auto_set_column_width([0])  # Adjust column width

    # Add a title to the table
    axs[1, 1].set_title('Results', fontsize=12, pad=20)  # Use 'pad' to adjust the distance from the top

    # Optionally add grid lines
    axs[1, 1].grid(False)

    # Optionally add a border around the table
    for key, cell in table.get_celld().items():
        cell.set_linewidth(0.5)

    plt.show()

    plt.tight_layout()
    plt.show()

    time_end = time.time()
    print("Time: ", time_end - time_start)