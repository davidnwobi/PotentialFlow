import numpy as np
import elementary_flows
from flow_field import FlowField
from multiprocessing import freeze_support
import random
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

    # x = np.linspace(X_NEG_LIMIT, X_POS_LIMIT, num=NO_OF_POINTS)
    # y = np.linspace(Y_NEG_LIMIT, Y_POS_LIMIT, num=NO_OF_POINTS)
    # #
    # velocity = 10
    # radius = 1
    # kappa = 2 * np.pi * velocity * radius ** 2 # Known solution for a cylinder
    # vortex_strength = 4 * np.pi * velocity * radius # Known solution for a cylinder
    #
    # alpha = np.pi/30 # Angle of attack
    # u1 = elementary_flows.UniformFlow(horizontal_vel=velocity * np.cos(alpha), vertical_vel=velocity * np.sin(alpha))
    # v1 = elementary_flows.Vortex(x_pos=0, y_pos=0, circulation=vortex_strength)
    # d1 = elementary_flows.Doublet(x_pos=0, y_pos=0, kappa=kappa)
    #
    # flow = FlowField([v1, u1, d1], **plotting_kwargs)
    # flow.plot_flow_from_stream_function(x, y).show()
    # flow.plot_velocity(x, y).show()
    #
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

    NO_OF_POINTS = 1000
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

    for i in range(500):
        flows.append(elementary_flows.Vortex(x_pos=random.uniform(-25, 25), y_pos=random.uniform(-25, 25), circulation=random.uniform(-10, 10)))


    flow = FlowField(flows, **plotting_kwargs)
    flow.plot_flow_from_stream_function(x, y).show()
    flow.plot_velocity(x, y).show()









