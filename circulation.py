import numpy as np
import scipy as sp
from tuple_collections import Ellipse, EllipseProperties, FlowFieldProperties


def compute_ellipse_and_circulation(flow_field: FlowFieldProperties, ellipse_def: Ellipse, divsions=100):
    t = np.linspace(0, 2 * np.pi, divsions)

    x_ellipse = ellipse_def.a * np.cos(t) + ellipse_def.x0
    y_ellipse = ellipse_def.b * np.sin(t) + ellipse_def.y0

    u_ellipse = sp.interpolate.RectBivariateSpline(flow_field.y, flow_field.x, flow_field.u).ev(y_ellipse, x_ellipse)
    v_ellipse = sp.interpolate.RectBivariateSpline(flow_field.y, flow_field.x, flow_field.v).ev(y_ellipse, x_ellipse)

    circulation = -(np.trapz(u_ellipse, x_ellipse) + np.trapz(v_ellipse, y_ellipse))
    return EllipseProperties(x_ellipse, y_ellipse, u_ellipse, v_ellipse, circulation)
