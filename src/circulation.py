import numpy as np
import scipy as sp
from . import data_collections as dc


def compute_ellipse_and_circulation(flow_field: dc.FlowFieldProperties, ellipse_def: dc.Ellipse, divsions: int=100):
    '''
    Computes the circulation around an ellipse placed in a flow field.

    Parameters
    ----------

    flow_field : FlowFieldProperties
    The flow field properties.

    ellipse_def : Ellipse
    The ellipse definition.

    divsions : int
    The number of divisions to use when discretizing the ellipse.

    Returns
    -------

    EllipseProperties
    The ellipse properties.

    '''

    t = np.linspace(0, 2 * np.pi, divsions)

    x_ellipse = ellipse_def.a * np.cos(t) + ellipse_def.x0
    y_ellipse = ellipse_def.b * np.sin(t) + ellipse_def.y0

    u_ellipse = sp.interpolate.RectBivariateSpline(flow_field.y, flow_field.x, flow_field.u).ev(y_ellipse, x_ellipse)
    v_ellipse = sp.interpolate.RectBivariateSpline(flow_field.y, flow_field.x, flow_field.v).ev(y_ellipse, x_ellipse)

    circulation = -(np.trapz(u_ellipse, x_ellipse) + np.trapz(v_ellipse, y_ellipse))
    return dc.EllipseProperties(x_ellipse, y_ellipse, u_ellipse, v_ellipse, circulation)
