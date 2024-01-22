# PANEL METHOD GEOMETRY
# Written by: JoshTheEngineer
# Started: 01/24/19
# Updated: 01/24/19 - Started code in MATLAB
#                   - Works as expected
#          02/03/19 - Transferred code to Python
#                   - Works as expected

import numpy as np
import matplotlib.pyplot as plt
from .tuple_collections import Geometry, PanelizedGeometry
from panel_generator import PanelGenerator

numB = 100  # Number of boundary points
tO = 22.5  # Angle offset [deg]
load = 'Circle'  # Load circle or airfoil
panelized_geometry = None
numPan = None
if load == 'Circle':  # If circle is selected
    theta = np.linspace(0, 360, numB)  # Angles to compute boundary points [deg]
    theta = theta + tO  # Add angle offset [deg]
    theta = theta * (np.pi / 180)  # Convert angle to radians [rad]

    # Boundary points
    XB = np.cos(theta)  # X value of boundary points
    YB = np.sin(theta)  # Y value of boundary points

    # Number of panels
    numPan = len(XB) - 1  # Number of panels

    circle = Geometry(XB, YB, 0)
    panelized_geometry = PanelGenerator.compute_geometric_quantities(circle)

else:  # If airfoil is selected
    # Load airfoil data
    data = np.loadtxt('naca2412.txt')
    XB = data[:, 0]
    YB = data[:, 1]

    # Number of panels
    numPan = len(XB) - 1
    airfoil = Geometry(XB, YB, 0)
    panelized_geometry = PanelGenerator.compute_geometric_quantities(airfoil)

# %% PLOTTING

# Dashed circle defined
T = np.linspace(0, 2 * np.pi, 1000)  # Angle array to compute circle
x = np.cos(T)  # Circle X points
y = np.sin(T)  # Circle Y points

# Plot the paneled geometry
fig = plt.figure(1)  # Create figure
plt.cla()  # Get ready for plotting
if (load == 'Circle'):  # If circle is selected
    plt.plot(x, y, 'k--')  # Plot actual circle outline
plt.fill(XB, YB, 'k')  # Plot polygon (circle or airfoil)
X = np.zeros(2)  # Initialize panel X variable
Y = np.zeros(2)  # Initialize panel Y variable
for i in range(numPan):  # Loop over all panels
    X[0] = panelized_geometry.control_points_x_cor[i]  # Panel starting X point
    X[1] = panelized_geometry.control_points_x_cor[i] + panelized_geometry.panel_length[i] * np.cos(
        panelized_geometry.panel_normal_angle[i])
    Y[0] = panelized_geometry.control_points_y_cor[i]  # Panel starting Y point
    Y[1] = panelized_geometry.control_points_y_cor[i] + panelized_geometry.panel_length[i] * np.sin(
        panelized_geometry.panel_normal_angle[i])
    if (i == 0):  # For first panel
        plt.plot(X, Y, 'b-', label='First Panel')  # Plot the first panel normal vector
    elif (i == 1):  # For second panel
        plt.plot(X, Y, 'g-', label='Second Panel')  # Plot the second panel normal vector
    else:  # For every other panel
        plt.plot(X, Y, 'r-')  # Plot the panel normal vector
plt.xlabel('X-Axis')  # Set X-label
plt.ylabel('Y-Axis')  # Set Y-label
plt.title('Panel Geometry')  # Set title
plt.axis('equal')  # Set axes equal
plt.legend()  # Plot legend
plt.show()  # Display plot
