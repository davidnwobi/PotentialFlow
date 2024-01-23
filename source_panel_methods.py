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


def COMPUTE_IJ_SPM(XC, YC, XB, YB, phi, S):
    # Number of panels
    numPan = len(XC)  # Number of panels/control points

    # Initialize arrays
    I = np.zeros([numPan, numPan])  # Initialize I integral matrix
    J = np.zeros([numPan, numPan])  # Initialize J integral matrix

    # Compute integral
    for i in range(numPan):  # Loop over i panels
        for j in range(numPan):  # Loop over j panels
            if (j != i):  # If the i and j panels are not the same
                # Compute intermediate values
                A = -(XC[i] - XB[j]) * np.cos(phi[j]) - (YC[i] - YB[j]) * np.sin(phi[j])  # A term
                B = (XC[i] - XB[j]) ** 2 + (YC[i] - YB[j]) ** 2  # B term
                Cn = np.sin(phi[i] - phi[j])  # C term (normal)
                Dn = -(XC[i] - XB[j]) * np.sin(phi[i]) + (YC[i] - YB[j]) * np.cos(phi[i])  # D term (normal)
                Ct = -np.cos(phi[i] - phi[j])  # C term (tangential)
                Dt = (XC[i] - XB[j]) * np.cos(phi[i]) + (YC[i] - YB[j]) * np.sin(phi[i])  # D term (tangential)
                E = np.sqrt(B - A ** 2)  # E term
                if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(
                        E)):  # If E term is 0 or complex or a NAN or an INF
                    I[i, j] = 0  # Set I value equal to zero
                    J[i, j] = 0  # Set J value equal to zero
                else:
                    # Compute I (needed for normal velocity), Ref [1]
                    term1 = 0.5 * Cn * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B)  # First term in I equation
                    term2 = ((Dn - A * Cn) / E) * (
                            math.atan2((S[j] + A), E) - math.atan2(A, E))  # Second term in I equation
                    I[i, j] = term1 + term2  # Compute I integral

                    # Compute J (needed for tangential velocity), Ref [2]
                    term1 = 0.5 * Ct * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B)  # First term in I equation
                    term2 = ((Dt - A * Ct) / E) * (
                            math.atan2((S[j] + A), E) - math.atan2(A, E))  # Second term in I equation
                    J[i, j] = term1 + term2  # Compute J integral

            # Zero out any problem values
            if (np.iscomplex(I[i, j]) or np.isnan(I[i, j]) or np.isinf(
                    I[i, j])):  # If I term is complex or a NAN or an INF
                I[i, j] = 0  # Set I value equal to zero
            if (np.iscomplex(J[i, j]) or np.isnan(J[i, j]) or np.isinf(
                    J[i, j])):  # If J term is complex or a NAN or an INF
                J[i, j] = 0  # Set J value equal to zero

    return I, J  # Return both I and J matrices


if __name__ == '__main__':
    freeze_support()
    numB = 9  # Number of boundary points
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
    AoA = 0
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
    Itest, Jtest = COMPUTE_IJ_SPM(panelized_geometry.xC, panelized_geometry.yC, X, Y, panelized_geometry.phi,
                                  panelized_geometry.S)
    print("I: ", I[0])
    print("Itest: ", Itest[0])
    gamma = compute_source_strengths(panelized_geometry, V, I)
    gamma_test = compute_source_strengths(panelized_geometry, V, Itest)
    sumLambda = sum(gamma * panelized_geometry.S)  # Check sum of source panel strengths
    print("Sum of L: ", sumLambda)
    print("Sum of L: ", sum(gamma_test * panelized_geometry.S))
    V_normal, V_tangential = compute_panel_velocities(panelized_geometry, gamma, V, I, J)

    Vx = V_tangential * np.cos(panelized_geometry.phi)
    Vy = V_tangential * np.sin(panelized_geometry.phi)
    # plt.figure(figsize=(12, 12), dpi=500)
    # plt.quiver(panelized_geometry.xC, panelized_geometry.yC, Vx, Vy)
    # plt.axis('equal')
    x, y = np.linspace(-2, 2, 100), np.linspace(-2, 2, 100)  # Create grid

    u, v = compute_grid_velocity(panelized_geometry, x, y, gamma, V, AoA)
    plt.figure(figsize=(12, 12), dpi=500)
    plt.streamplot(x, y, u, v, density=4, color='b')
    # plt.quiver(x, y, Vx, Vy)
    plt.xlim(-1, 2)
    plt.ylim(-1, 2)
    plt.fill(geometry.x, geometry.y, 'k')  # Plot polygon (circle or airfoil)
    plt.axis('equal')
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
    # plt.axis('equal')
    # plt.show()
