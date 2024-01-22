import numpy as np
from src import elementary_flows
from src.flow_field import FlowField
from multiprocessing import freeze_support
import matplotlib.pyplot as plt
import random
from src.data_collections import FlowFieldProperties, Ellipse, EllipseProperties
from src.circulation import compute_ellipse_and_circulation

if __name__ == '__main__':
    freeze_support()
    NO_OF_POINTS = 1000
    X_POS_LIMIT = 5
    Y_POS_LIMIT = 5
    X_NEG_LIMIT = -5
    Y_NEG_LIMIT = -5

    plotting_kwargs = {
        'X_NEG_LIMIT': X_NEG_LIMIT,
        'X_POS_LIMIT': X_POS_LIMIT,
        'Y_NEG_LIMIT': Y_NEG_LIMIT,
        'Y_POS_LIMIT': Y_POS_LIMIT,
        'STREAMLINE_DENSITY': 3,
        'STREAMLINE_COLOR': 'b',
        'CONTOR_LEVELS': 100,
        'CONTOR_COLOR': 'k',
        'FIGURE_SIZE': (12, 12),
        'DPI': 100,
        "CONTOUR_LABELS": True
    }

    x = np.linspace(X_NEG_LIMIT, X_POS_LIMIT, num=NO_OF_POINTS)
    y = np.linspace(Y_NEG_LIMIT, Y_POS_LIMIT, num=NO_OF_POINTS)
    #
    velocity = 10
    radius = 1
    kappa = 2 * np.pi * velocity * radius ** 2  # Known solution for a cylinder
    vortex_strength = 4 * np.pi * velocity * radius  # Known solution for a cylinder

    alpha =0* np.pi / 30  # Angle of attack
    u1 = elementary_flows.UniformFlow(horizontal_vel=velocity * np.cos(alpha), vertical_vel=velocity * np.sin(alpha))
    v1 = elementary_flows.Vortex(x_pos=0, y_pos=0, circulation=vortex_strength)
    d1 = elementary_flows.Doublet(x_pos=0, y_pos=0, kappa=kappa)
    s1= elementary_flows.Source(x_pos=0, y_pos=0, strength=10)
    flow = FlowField([s1], **plotting_kwargs)
    X, Y = np.meshgrid(x, y)
    stream_function = np.array(flow.stream_function(X, Y))
    velocity_field = flow.velocity(X, Y)
    fig = flow.plot_flow_from_stream_function(x, y)

    # Create ellipse and compute circulation
    flow_properties = FlowFieldProperties(x, y, velocity_field[0], velocity_field[1])
    ellipse = Ellipse(x0=-2, y0=-2, a=1, b=1)
    ellipse_properties = compute_ellipse_and_circulation(flow_properties, ellipse, divsions=1000)

    # Plot ellipse and circulation
    plt.plot(ellipse_properties.x_cor, ellipse_properties.y_cor, color='r', linewidth=5)
    plt.quiver(ellipse_properties.x_cor, ellipse_properties.y_cor, ellipse_properties.u, ellipse_properties.v,
               color='r', scale=10000)
    plt.show()
    print(f"circulation: {ellipse_properties.circulation}")
