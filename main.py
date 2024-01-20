import numpy as np
import elementary_flows
from flow_field import FlowField


NO_OF_POINTS = 1000
X_POS_LIMIT = 5
Y_POS_LIMIT = 5
X_NEG_LIMIT = -5
Y_NEG_LIMIT = -5

x = np.linspace(X_NEG_LIMIT, X_POS_LIMIT, num=NO_OF_POINTS)
y = np.linspace(Y_NEG_LIMIT, Y_POS_LIMIT, num=NO_OF_POINTS)
velocity = 10
radius = 1
kappa = 2 * np.pi * velocity * radius ** 2
vortex_strength = 4 * np.pi * velocity * radius*0.5

alpha = np.pi/30
u1 = elementary_flows.UniformFlow(horizontal_vel=velocity * np.cos(alpha), vertical_vel=velocity * np.sin(alpha))
v1 = elementary_flows.Vortex(x_pos=0, y_pos=0, circulation=vortex_strength)
d1 = elementary_flows.Doublet(x_pos=0, y_pos=0, kappa=kappa)


plotting_kwargs = {
    'X_NEG_LIMIT': X_NEG_LIMIT,
    'X_POS_LIMIT': X_POS_LIMIT,
    'Y_NEG_LIMIT': Y_NEG_LIMIT,
    'Y_POS_LIMIT': Y_POS_LIMIT,
    'STREAMLINE_DENSITY': 3,
    'STREAMLINE_COLOR': 'b',
    'CONTOR_LEVELS': 500,
    'CONTOR_COLOR': 'k',
    'FIGURE_SIZE': (12, 12),
    'DPI': 100,
    "CONTOUR_LABELS": True
}
flow = FlowField([v1, u1, d1], **plotting_kwargs)
flow.plot(x, y).show()
