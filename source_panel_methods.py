import matplotlib.pyplot as plt
from src.panel_generator import PanelGenerator
from src.geometric_integrals import compute_panel_geometric_integrals, compute_grid_geometric_integrals
import src.data_collections as dc
from src.airfoil_generator import generate_four_digit_NACA
from src.elementary_flows import *
from src.flow_field import *
from src.source_panel_methods_funcs import *
from src.plotting import plot_panelized_geometry
from multiprocessing import freeze_support
import numpy as np
import math

# FUNCTION - COMPUTE I AND J GEOMETRIC INTEGRALS FOR SOURCE PANEL METHOD
# Written by: JoshTheEngineer
# YouTube   : www.youtube.com/joshtheengineer
# Website   : www.joshtheengineer.com
# Started   : 02/03/19 - Transferred from MATLAB to Python
#                      - Works as expected
#           : 04/28/20 - Fixed E value error handling
#
# PURPOSE
# - Compute the integral expression for constant strength source panels
# - Source panel strengths are constant, but can change from panel to panel
# - Geometric integral for panel-normal    : I(ij)
# - Geometric integral for panel-tangential: J(ij)
#
# REFERENCES
# - [1]: Normal Geometric Integral SPM, I(ij)
#           Link: https://www.youtube.com/watch?v=76vPudNET6U
# - [2]: Tangential Geometric Integral SPM, J(ij)
#           Link: https://www.youtube.com/watch?v=JRHnOsueic8
#
# INPUTS
# - XC  : X-coordinate of control points
# - YC  : Y-coordinate of control points
# - XB  : X-coordinate of boundary points
# - YB  : Y-coordinate of boundary points
# - phi : Angle between positive X-axis and interior of panel
# - S   : Length of panel
#
# OUTPUTS
# - I   : Value of panel-normal integral (Eq. 3.163 in Anderson or Ref [1])
# - J   : Value of panel-tangential integral (Eq. 3.165 in Anderson or Ref [2])

import numpy as np
import math as math




if __name__ == '__main__':
    freeze_support()
    numB = 100  # Number of boundary points
    tO = 0  # Angle offset [deg]
    load = 'CircleS'  # Load circle or airfoil
    geometry = None
    panelized_geometry = None
    numPan = None

    theta = np.linspace(0, 360, numB)  # Angles to compute boundary points [deg]
    theta = theta + tO  # Add angle offset [deg]
    theta = theta * (np.pi / 180)  # Convert angle to radians [rad]

    # Boundary points
    XB = np.cos(theta)  # X value of boundary points
    YB = np.sin(theta)  # Y value of boundary points

    # Number of panels
    numPan = len(XB) - 1  # Number of panels
    AoA = 10
    geometry = dc.Geometry(XB, YB, AoA)
    if not load == 'Circle':  # If circle is selected
        XB, YB = generate_four_digit_NACA("0012", numB, 1)
        # Number of panels
        geometry = dc.Geometry(XB, YB, AoA)
    panelized_geometry = PanelGenerator.compute_geometric_quantities(geometry)
    # plot_panelized_geometry(geometry, panelized_geometry)

    V = 10
    X = panelized_geometry.xC - panelized_geometry.S / 2 * np.cos(panelized_geometry.phi)
    Y = panelized_geometry.yC - panelized_geometry.S / 2 * np.sin(panelized_geometry.phi)
    I, J = compute_panel_geometric_integrals(panelized_geometry)

    gamma = compute_source_strengths(panelized_geometry, V, I)

    sumLambda = sum(gamma * panelized_geometry.S)  # Check sum of source panel strengths
    print("Sum of L: ", sumLambda)
    V_normal, V_tangential = compute_panel_velocities(panelized_geometry, gamma, V, I, J)
    # plt.figure(figsize=(12, 12), dpi=500)
    # plt.quiver(panelized_geometry.xC, panelized_geometry.yC, Vx, Vy)
    # plt.axis('equal')
    x, y = np.linspace(-0.5, 1.5, 200), np.linspace(-0.5, 0.5, 200)  # Create grid

    u, v = compute_grid_velocity(panelized_geometry, x, y, gamma, V, AoA)
    plt.figure(figsize=(12, 12), dpi=500)
    plt.streamplot(x, y, u, v, density=2, color='b')
    # plt.quiver(x, y, Vx, Vy)
    plt.xlim(-0.5, 1.5)
    plt.ylim(-0.5, 0.5)
    plt.fill(geometry.x, geometry.y, 'k')  # Plot polygon (circle or airfoil)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()

    # plot panel geometry

    plt.figure(figsize=(12, 12), dpi=500)

    X = np.append(X, X[0])
    Y = np.append(Y, Y[0])
    panel_vector_X = panelized_geometry.xC + panelized_geometry.S / 2 * np.cos(panelized_geometry.delta)
    panel_vector_Y = panelized_geometry.yC + panelized_geometry.S / 2 * np.sin(panelized_geometry.delta)
    # for i in range(len(panelized_geometry.S)):
    #     plt.plot([X[i], X[i + 1]], [Y[i], Y[i + 1]], 'k')
    #     plt.plot([X[i]], [Y[i]], 'ro')
    #     plt.plot([panelized_geometry.xC[i]], [panelized_geometry.yC[i]], 'bo')
    #     plt.text(panelized_geometry.xC[i], panelized_geometry.yC[i], str(i))
    #     plt.plot([panelized_geometry.xC[i], panel_vector_X[i]], [panelized_geometry.yC[i], panel_vector_Y[i]], 'k')
    #     plt.plot([panel_vector_X[i]], [panel_vector_Y[i]], 'go')


    # plt.show()
