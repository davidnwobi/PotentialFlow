import numpy as np
from src import elementary_flows
from flow_field import FlowField
from multiprocessing import freeze_support
import matplotlib.pyplot as plt
import random
from tuple_collections import FlowFieldProperties, Ellipse, EllipseProperties
from circulation import compute_ellipse_and_circulation

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

    flow = FlowField([v1, d1, u1], **plotting_kwargs)
    flow.plot_flow_from_stream_function(x, y)
    flow.plot_velocity(x, y)


    # Rankine Oval
    plotting_kwargs2 = {
        'CONTOR_LEVELS': 50,
    }
    v1 = elementary_flows.Source(x_pos=-2, y_pos=0, strength=10)
    v2 = elementary_flows.Source(x_pos=2, y_pos=0, strength=-10)# Negative strength is a sink
    u1 = elementary_flows.UniformFlow(horizontal_vel=1, vertical_vel=0)

    flow = FlowField([v1, v2, u1], **plotting_kwargs2)
    flow.plot_flow_from_stream_function(x, y).show()
    flow.plot_velocity(x, y).show()



    # Kelvin's Oval

    plotting_kwargs2 = {
        'CONTOR_LEVELS': 50,
    }
    v1 = elementary_flows.Vortex(x_pos=0, y_pos=2, circulation=10)
    v2 = elementary_flows.Vortex(x_pos=0, y_pos=-2, circulation=-10)
    u1 = elementary_flows.UniformFlow(horizontal_vel=1, vertical_vel=0)
    flow = FlowField([v1, v2, u1], **plotting_kwargs2)
    flow.plot_flow_from_stream_function(x, y).show()
    flow.plot_velocity(x, y).show()


    # Random flow field to test multiprocessing
    NO_OF_POINTS = 50
    X_POS_LIMIT = 30
    Y_POS_LIMIT = 30
    X_NEG_LIMIT = -30
    Y_NEG_LIMIT = -30

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
    u1 = elementary_flows.UniformFlow(horizontal_vel=1, vertical_vel=0)
    flow_types = [elementary_flows.Source, elementary_flows.Vortex, elementary_flows.Doublet]
    flows = []
    flows.append(u1)

    #  Define random elementary flows
    for i in range(50):
        flow_type = random.choice(flow_types)
        if flow_type == elementary_flows.Source:
            flows.append(flow_type(x_pos=random.uniform(-20, 20), y_pos=random.uniform(-20, 20),
                                   strength=random.uniform(-10, 10)))
        elif flow_type == elementary_flows.Vortex:
            flows.append(flow_type(x_pos=random.uniform(-20, 20), y_pos=random.uniform(-20, 20),
                                   circulation=random.uniform(-10, 10)))
        elif flow_type == elementary_flows.Doublet:
            flows.append(
                flow_type(x_pos=random.uniform(-20, 20), y_pos=random.uniform(-20, 20), kappa=random.uniform(-10, 10)))

    # Create random flow field
    # flow = FlowField(flows, **plotting_kwargs)
    # X, Y = np.meshgrid(x, y)
    # velocity_field = flow.velocity(X, Y)
    #
    # # Plot velocity field
    # fig = plt.figure(figsize=plotting_kwargs.get("FIGURE_SIZE", (12, 12)), dpi=plotting_kwargs.get("DPI", 100))
    # plt.quiver(X, Y, velocity_field[0], velocity_field[1])
