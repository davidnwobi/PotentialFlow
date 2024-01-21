import numpy as np
import elementary_flows
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
    X, Y = np.meshgrid(x, y)
    velocity_field = flow.velocity(X, Y)
    fig = flow.plot_velocity(x, y)

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

    # # Rankine Oval
    # plotting_kwargs2 = {
    #     'CONTOR_LEVELS': 50,
    # }
    # v1 = elementary_flows.Source(x_pos=-2, y_pos=0, strength=10)
    # v2 = elementary_flows.Source(x_pos=2, y_pos=0, strength=-10)# Negative strength is a sink
    # u1 = elementary_flows.UniformFlow(horizontal_vel=1, vertical_vel=0)
    #
    # flow = FlowField([v1, v2, u1], **plotting_kwargs2)
    # flow.plot_flow_from_stream_function(x, y).show()
    # flow.plot_velocity(x, y).show()
    #
    #
    #
    # # Kelvin's Oval
    #
    # plotting_kwargs2 = {
    #     'CONTOR_LEVELS': 50,
    # }
    # v1 = elementary_flows.Vortex(x_pos=0, y_pos=2, circulation=10)
    # v2 = elementary_flows.Vortex(x_pos=0, y_pos=-2, circulation=-10)
    # u1 = elementary_flows.UniformFlow(horizontal_vel=1, vertical_vel=0)
    # flow = FlowField([v1, v2, u1], **plotting_kwargs2)
    # flow.plot_flow_from_stream_function(x, y).show()
    # flow.plot_velocity(x, y).show()

    # NO_OF_POINTS = 50
    # X_POS_LIMIT = 30
    # Y_POS_LIMIT = 30
    # X_NEG_LIMIT = -30
    # Y_NEG_LIMIT = -30
    #
    # plotting_kwargs = {
    #     'X_NEG_LIMIT': X_NEG_LIMIT,
    #     'X_POS_LIMIT': X_POS_LIMIT,
    #     'Y_NEG_LIMIT': Y_NEG_LIMIT,
    #     'Y_POS_LIMIT': Y_POS_LIMIT,
    #     'STREAMLINE_DENSITY': 3,
    #     'STREAMLINE_COLOR': 'b',
    #     'CONTOR_LEVELS': 100,
    #     'CONTOR_COLOR': 'k',
    #     'FIGURE_SIZE': (12, 12),
    #     'DPI': 100,
    #     "CONTOUR_LABELS": True
    # }
    #
    # x = np.linspace(X_NEG_LIMIT, X_POS_LIMIT, num=NO_OF_POINTS)
    # y = np.linspace(Y_NEG_LIMIT, Y_POS_LIMIT, num=NO_OF_POINTS)
    # u1 = elementary_flows.UniformFlow(horizontal_vel=1, vertical_vel=0)
    # flow_types = [elementary_flows.Source, elementary_flows.Vortex, elementary_flows.Doublet]
    # flows = []
    # flows.append(u1)
    #
    # #  Define random elementary flows
    # for i in range(50):
    #     flow_type = random.choice(flow_types)
    #     if flow_type == elementary_flows.Source:
    #         flows.append(flow_type(x_pos=random.uniform(-20, 20), y_pos=random.uniform(-20, 20),
    #                                strength=random.uniform(-10, 10)))
    #     elif flow_type == elementary_flows.Vortex:
    #         flows.append(flow_type(x_pos=random.uniform(-20, 20), y_pos=random.uniform(-20, 20),
    #                                circulation=random.uniform(-10, 10)))
    #     elif flow_type == elementary_flows.Doublet:
    #         flows.append(
    #             flow_type(x_pos=random.uniform(-20, 20), y_pos=random.uniform(-20, 20), kappa=random.uniform(-10, 10)))
    #
    # # Create random flow field
    # flow = FlowField(flows, **plotting_kwargs)
    # X, Y = np.meshgrid(x, y)
    # velocity_field = flow.velocity(X, Y)
    #
    # # Plot velocity field
    # fig = plt.figure(figsize=plotting_kwargs.get("FIGURE_SIZE", (12, 12)), dpi=plotting_kwargs.get("DPI", 100))
    # plt.quiver(X, Y, velocity_field[0], velocity_field[1])
    #
    # # Create ellipse and compute circulation
    # flow_properties = FlowFieldProperties(x, y, velocity_field[0], velocity_field[1])
    # ellipse = Ellipse(x0=0, y0=0, a=10, b=5)
    # ellipse_properties = compute_ellipse_and_circulation(flow_properties, ellipse, divsions=100)
    #
    # # Plot ellipse and circulation
    # plt.plot(ellipse_properties.x_cor, ellipse_properties.y_cor)
    # plt.quiver(ellipse_properties.x_cor, ellipse_properties.y_cor, ellipse_properties.u, ellipse_properties.v, color='r', scale=50)
    # plt.show()
    # print(f"circulation: {ellipse_properties.circulation}")
