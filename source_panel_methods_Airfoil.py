import matplotlib.pyplot as plt
from src.panel_generator import PanelGenerator
from src.geometric_integrals import compute_panel_geometric_integrals
import src.data_collections as dc
from src.airfoil_generator import generate_four_digit_NACA
from src.source_panel_methods_funcs import *
from aeropy import xfoil_module as xf

import numpy as np
import pandas as pd

if __name__ == '__main__':
    airfoil = '2412'
    AoA = 8
    res = xf.find_pressure_coefficients(airfoil='naca' + airfoil, alpha=AoA, NACA=True, )

    xfoil_cp = pd.DataFrame(res)
    xfoil_cp_upp = xfoil_cp[xfoil_cp['y'] >= 0]
    xfoil_cp_low = xfoil_cp[xfoil_cp['y'] <= 0]

    numB = 25  # Number of boundary points
    num_grid = 25
    X_NEG_LIMIT = -0.5
    X_POS_LIMIT = 1.5
    Y_NEG_LIMIT = -1
    Y_POS_LIMIT = 1

    load = 'Circle2'  # Load circle or airfoil

    XB, YB = generate_four_digit_NACA(airfoil, numB, 1)

    V = 1
    geometry = dc.Geometry(XB, YB, AoA)
    panelized_geometry = PanelGenerator.compute_geometric_quantities(geometry)

    I, J = compute_panel_geometric_integrals(panelized_geometry)
    gamma = compute_source_strengths(panelized_geometry, V, I)
    sumLambda = sum(gamma * panelized_geometry.S)  # Check sum of source panel strengths
    print("Sum of L: ", sumLambda)

    V_normal, V_tangential = compute_panel_velocities(panelized_geometry, gamma, V, I, J)
    x, y = np.linspace(X_NEG_LIMIT, X_POS_LIMIT, num_grid), np.linspace(Y_NEG_LIMIT, Y_POS_LIMIT, num_grid)  # Create grid
    u, v = compute_grid_velocity(panelized_geometry, x, y, gamma, V, AoA)

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
    spm_CP = pd.DataFrame({
        'x': panelized_geometry.xC,
        'y': panelized_geometry.yC,
        'Cp': Cp
    })
    spm_CP_upp = spm_CP[spm_CP['y'] >= 0]
    spm_CP_low = spm_CP[spm_CP['y'] <= 0]

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
        if i == 0:
            axs[0, 0].plot([panelized_geometry.xC[i], panel_vector_X[i]], [panelized_geometry.yC[i], panel_vector_Y[i]],
                           'k', label='First Panel')
        if i == 1:
            axs[0, 0].plot([panelized_geometry.xC[i], panel_vector_X[i]], [panelized_geometry.yC[i], panel_vector_Y[i]],
                           'k', label='Second Panel')
        else:
            axs[0, 0].plot([panelized_geometry.xC[i], panel_vector_X[i]], [panelized_geometry.yC[i], panel_vector_Y[i]],
                           'k')
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
    axs[1, 0].plot(xfoil_cp_upp['x'], xfoil_cp_upp['Cp'], 'r', label='Upper Surface Xfoil')
    axs[1, 0].plot(xfoil_cp_low['x'], xfoil_cp_low['Cp'], 'b', label='Lower Surface Xfoil')
    axs[1, 0].plot(spm_CP_upp['x'], spm_CP_upp['Cp'], 'ro', label='Upper Surface SPM')
    axs[1, 0].plot(spm_CP_low['x'], spm_CP_low['Cp'], 'bo', label='Lower Surface SPM')
    axs[1, 0].legend()
    axs[1, 0].invert_yaxis()
    axs[1, 0].set_xlabel('x/c')
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
