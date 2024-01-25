import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as colors
from matplotlib import ticker
from src.panel_generator import PanelGenerator
from src.circulation import compute_ellipse_and_circulation
from src.vortex_panel_methods_funcs import *
from src.airfoil_generator import generate_four_digit_NACA
from aeropy import xfoil_module as xf
import numpy as np
import pandas as pd
import time

if __name__ == '__main__':
    airfoil = '2412'
    AoA = 6
    res = None
    x_foil_cl = 0

    try:
        res = xf.find_pressure_coefficients(airfoil='naca' + airfoil, alpha=AoA, NACA=True, )
        x_foil_cl = xf.find_coefficients(airfoil='naca' + airfoil, alpha=AoA, NACA=True, )['CL']
        xfoil_cp = pd.DataFrame(res)
        xfoil_cp_upp = xfoil_cp[xfoil_cp['y'] >= 0]
        xfoil_cp_low = xfoil_cp[xfoil_cp['y'] <= 0]
    except:
        print('XFOIL Error')
        res = xf.find_pressure_coefficients(airfoil='naca' + airfoil, alpha=0, NACA=True, )
        xfoil_cp = pd.DataFrame(res)
        xfoil_cp_upp = xfoil_cp[xfoil_cp['y'] >= 0]
        xfoil_cp_low = xfoil_cp[xfoil_cp['y'] <= 0]

    num_grid = 50  # Change this if it is too slow or you run out of memory
    X_NEG_LIMIT = -0.5
    X_POS_LIMIT = 1.25
    Y_NEG_LIMIT = -0.5
    Y_POS_LIMIT = 0.5

    load = 'Circle2'  # Load circle or airfoil

    time_start = time.time()
    # %% THE ABOVE CODE DOES NOT PROVIDE THE RIGHT BOUNDARY POINTS FOR THE AIRFOIL. LOAD THE AIRFOIL DATA FROM A FILE

    # %% PANEL METHOD GEOMETRY SETUP
    XB, YB = generate_four_digit_NACA(num_NACA=airfoil, num_points=30, chord_length=1)
    numB = len(XB)
    V = 1
    geometry = dc.Geometry(XB, YB, AoA)
    panelized_geometry = PanelGenerator.compute_geometric_quantities(geometry)
    x, y = np.linspace(X_NEG_LIMIT, X_POS_LIMIT, num_grid), np.linspace(Y_NEG_LIMIT, Y_POS_LIMIT,
                                                                        num_grid)  # Create grid

    X = geometry.x
    Y = geometry.y
    panel_normal_vector_X = panelized_geometry.xC + panelized_geometry.S / 2 * np.cos(panelized_geometry.delta)
    panel_normal_vector_Y = panelized_geometry.yC + panelized_geometry.S / 2 * np.sin(panelized_geometry.delta)

    # %% VORTEX PANEL METHOD
    V_normal, V_tangential, gamma, u, v = run_vortex_panel_method(panelized_geometry, V, AoA, x, y)
    sumGamma = np.sum(gamma * panelized_geometry.S)
    print('Sum of Vortex Circulation: ', sumGamma)

    # %% Pressure Coefficient for the grid
    local_v = np.sqrt(u ** 2 + v ** 2)
    cp = 1 - (local_v / V) ** 2

    # %% Pressure Coefficient for the panels
    panel_velocities = V_tangential ** 2
    Cp = 1 - panel_velocities / V ** 2
    vpm_CP = pd.DataFrame({
        'x': panelized_geometry.xC,
        'y': panelized_geometry.yC,
        'Cp': Cp,
        'V': panel_velocities
    })
    vpm_CP_upp = vpm_CP[vpm_CP['y'] >= 0]
    vpm_CP_low = vpm_CP[vpm_CP['y'] <= 0]

    # %% Lift and Drag Coefficients
    CN = -Cp * np.sin(panelized_geometry.beta) * panelized_geometry.S  # Normal coefficient
    CA = -Cp * np.cos(panelized_geometry.beta) * panelized_geometry.S  # Axial coefficient

    CL = np.sum(CN * np.cos(AoA * np.pi / 180)) - np.sum(CA * np.sin(AoA * np.pi / 180))  # Lift coefficient
    CD = np.sum(CN * np.sin(AoA * np.pi / 180)) + np.sum(CA * np.cos(AoA * np.pi / 180))  # Drag coefficient

    airfoil_ellipse = dc.Ellipse(0.5, 0, 0.6, 0.2)
    flowfieldproperties = dc.FlowFieldProperties(x=x, y=y, u=u, v=v)
    circulation = compute_ellipse_and_circulation(flowfieldproperties, airfoil_ellipse, 1000)

    time_end = time.time()
    print('Time elapsed: ', time_end - time_start)

    # %% PLOTTING
    fig, axs = plt.subplots(2, 3, figsize=(25, 12), dpi=100)
    fig.suptitle(f'Vortex Panel Method Results for NACA {airfoil} at {AoA} degrees AoA', fontsize=16)
    # Panel Geometry
    axs[0, 0].set_title('Panel Geometry')
    for i in range(len(panelized_geometry.S)):
        axs[0, 0].plot([X[i], X[i + 1]], [Y[i], Y[i + 1]], 'k')
        axs[0, 0].plot([X[i]], [Y[i]], 'ro', markersize=0.5)
        axs[0, 0].plot([panelized_geometry.xC[i]], [panelized_geometry.yC[i]], 'bo', markersize=0.5)
        if i == 0:
            axs[0, 0].plot([panelized_geometry.xC[i], panel_normal_vector_X[i]],
                           [panelized_geometry.yC[i], panel_normal_vector_Y[i]],
                           'k', label='First Panel')
        if i == 1:
            axs[0, 0].plot([panelized_geometry.xC[i], panel_normal_vector_X[i]],
                           [panelized_geometry.yC[i], panel_normal_vector_Y[i]],
                           'k', label='Second Panel')
        else:
            axs[0, 0].plot([panelized_geometry.xC[i], panel_normal_vector_X[i]],
                           [panelized_geometry.yC[i], panel_normal_vector_Y[i]],
                           'k')
        axs[0, 0].plot([panel_normal_vector_X[i]], [panel_normal_vector_Y[i]], 'go', markersize=0.5)
    axs[0, 0].set_xlabel('X')
    axs[0, 0].set_ylabel('Y')
    axs[0, 0].set_aspect('equal', adjustable='box')

    # Streamlines
    axs[0, 1].set_title('Streamlines')
    axs[0, 1].streamplot(x, y, u, v, density=1.5, color='b')
    axs[0, 1].set_xlim(X_NEG_LIMIT, X_POS_LIMIT)
    axs[0, 1].set_ylim(Y_NEG_LIMIT, Y_POS_LIMIT)
    axs[0, 1].fill(geometry.x, geometry.y, 'k')  # Plot polygon (circle or airfoil)
    # draw_circle at the trailing edge
    axs[0, 1].add_patch(
        patches.Circle(
            (geometry.x[0], geometry.y[0]),
            0.1,
            fill=False,
            color='r',
            linewidth=3
        ))
    axs[0, 1].annotate('Kutta Condition', xy=(geometry.x[0], geometry.y[0]), xycoords='data',
                       xytext=(.99, -.1), textcoords='axes fraction',
                       va='top', ha='left',
                       arrowprops=dict(facecolor='black', shrink=0.05))
    axs[0, 1].set_xlabel('X')
    axs[0, 1].set_ylabel('Y')
    axs[0, 1].set_aspect('equal', adjustable='box')

    # Pressure Contours
    axs[0, 2].set_title('Pressure Contours')
    x_m, y_m = np.meshgrid(x, y)
    cp_plot = axs[0, 2].pcolor(x_m, y_m, cp, cmap='jet', norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,
                                                                                vmin=cp.min(), vmax=cp.max()))
    ticks = -np.flip(np.geomspace(0.05, np.abs(cp.min()), 7))
    ticks = np.append(ticks, np.geomspace(0.01, cp.max(), 3))
    cbar = fig.colorbar(cp_plot, ax=axs[0, 2], label="Pressure Coefficient", ticks=ticks)
    cbar.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))

    axs[0, 2].fill(geometry.x, geometry.y, 'k')

    axs[0, 2].set_ylim(Y_NEG_LIMIT, Y_POS_LIMIT)
    axs[0, 2].set_xlim(X_NEG_LIMIT, X_POS_LIMIT)
    axs[0, 2].set_aspect('equal', adjustable='box')
    axs[0, 2].set_xlabel('X')
    axs[0, 2].set_ylabel('Y')
    axs[0, 2].set_aspect('equal', adjustable='box')

    # Cp Distribution
    axs[1, 0].set_title('Cp Distribution')
    axs[1, 0].invert_yaxis()
    axs[1, 0].plot(xfoil_cp_upp['x'], xfoil_cp_upp['Cp'], 'r', label='Upper Surface Xfoil')
    axs[1, 0].plot(xfoil_cp_low['x'], xfoil_cp_low['Cp'], 'b', label='Lower Surface Xfoil')
    axs[1, 0].plot(vpm_CP_upp['x'], vpm_CP_upp['Cp'], 'bo', label='Upper Surface VPM')
    axs[1, 0].plot(vpm_CP_low['x'], vpm_CP_low['Cp'], 'ro', label='Lower Surface VPM')
    axs[1, 0].legend()
    axs[1, 0].set_xlabel('x/c')
    axs[1, 0].set_ylabel('Cp')
    axs[1, 0].legend()

    # Pressure Vectors
    axs[1, 1].set_title('Pressure Vectors')
    dCPx = np.cos(panelized_geometry.beta) * Cp
    dCPy = np.sin(panelized_geometry.beta) * Cp
    positive_top = np.where(Cp >= 0, True, False) & np.where(panelized_geometry.yC >= 0, True, False)
    positive_bottom = np.where(Cp >= 0, True, False) & np.where(panelized_geometry.yC < 0, True, False)
    negative_top = np.where(Cp < 0, True, False) & np.where(panelized_geometry.yC >= 0, True, False)
    negative_bottom = np.where(Cp < 0, True, False) & np.where(panelized_geometry.yC < 0, True, False)
    axs[1, 1].fill(geometry.x, geometry.y, 'k')
    axs[1, 1].quiver(panelized_geometry.xC[positive_top], panelized_geometry.yC[positive_top], -dCPx[positive_top],
                     -dCPy[positive_top], pivot='tip', scale=10, color='r', width=0.0025, headwidth=2, headlength=4)
    axs[1, 1].quiver(panelized_geometry.xC[positive_bottom], panelized_geometry.yC[positive_bottom],
                     dCPx[positive_bottom], dCPy[positive_bottom], pivot='tail', scale=10, color='r', width=0.0025,
                     headwidth=2, headlength=4)
    axs[1, 1].quiver(panelized_geometry.xC[negative_top], panelized_geometry.yC[negative_top], -dCPx[negative_top],
                     -dCPy[negative_top], pivot='tail', scale=10, color='b', width=0.0025, headwidth=2, headlength=4)
    axs[1, 1].quiver(panelized_geometry.xC[negative_bottom], panelized_geometry.yC[negative_bottom],
                     dCPx[negative_bottom], dCPy[negative_bottom], pivot='tip', scale=10, color='b', width=0.0025,
                     headwidth=2, headlength=4)
    axs[1, 1].set_xlim(X_NEG_LIMIT, X_POS_LIMIT)
    axs[1, 1].set_ylim(Y_NEG_LIMIT, Y_POS_LIMIT)
    axs[1, 1].set_xlabel('x/c')
    axs[1, 1].set_ylabel('Cp')
    axs[1, 1].set_aspect('equal', adjustable='box')

    # Results Table
    axs[1, 2].axis('off')
    table_data = pd.DataFrame({
        "Sum of Vortex Circulation": [round(sumGamma, 6)],
        "cl (Calculated from CP)": [round(CL, 6)],
        'Circulation (Evaluated from grid velocities)': [round(circulation.circulation, 6)],
        'cl (Kutta-Joukowski)': [round(2 * circulation.circulation, 6)],
        'cl (XFOIL)': [round(x_foil_cl, 6)],

    })
    table_data = table_data.round(6)
    table_data = table_data.T
    row_labels = table_data.index.values
    # Customize the table appearance
    table = axs[1, 2].table(cellText=table_data.values,
                            rowLabels=row_labels,
                            loc='center',
                            cellLoc='center',
                            bbox=[0.5, 0.6, 0.18, 0.2])

    # Adjust font size and style
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.auto_set_column_width([0])  # Adjust column width

    # Add a title to the table
    axs[1, 2].set_title('Results', fontsize=12, pad=20)  # Use 'pad' to adjust the distance from the top

    # Optionally add grid lines
    axs[1, 2].grid(False)

    # Optionally add a border around the table
    for key, cell in table.get_celld().items():
        cell.set_linewidth(0.5)

    plt.show()
