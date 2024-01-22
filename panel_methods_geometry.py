# PANEL METHOD GEOMETRY
# Written by: JoshTheEngineer
# Started: 01/24/19
# Updated: 01/24/19 - Started code in MATLAB
#                   - Works as expected
#          02/03/19 - Transferred code to Python
#                   - Works as expected

import numpy as np
import matplotlib.pyplot as plt
from src.data_collections import Geometry, PanelizedGeometry
from src.panel_generator import PanelGenerator
from src.plotting import plot_panelized_geometry

numB = 10  # Number of boundary points
tO = 22.5  # Angle offset [deg]
load = 'Circles'  # Load circle or airfoil
geometry = None
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

    geometry = Geometry(XB, YB, 0)
    panelized_geometry = PanelGenerator.compute_geometric_quantities(geometry)

else:  # If airfoil is selected
    # Load airfoil data
    data = np.loadtxt('naca2412.txt')
    XB = data[:, 0]
    YB = data[:, 1]

    # Number of panels
    numPan = len(XB) - 1
    geometry = Geometry(XB, YB, 0)
    panelized_geometry = PanelGenerator.compute_geometric_quantities(geometry)

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


fig = plot_panelized_geometry(geometry, panelized_geometry)
if (load == 'Circle'):  # If circle is selected
    plt.plot(x, y, 'k--')  # Plot actual circle outline
plt.show()