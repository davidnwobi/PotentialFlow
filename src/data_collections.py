from collections import namedtuple
import numpy as np
from dataclasses import dataclass
import typing as tp


class FlowFieldProperties(namedtuple("FlowField", ['x', 'y', 'u', 'v'])):
    '''
    A class that represents a flow field properties.

    Attributes
    ----------

    x : np.ndarray
    The x-array coordinates of the flow field.

    y : np.ndarray
    The y-array coordinates of the flow field.

    u : np.ndarray
    The x-component of the velocity field.

    v : np.ndarray
    The y-component of the velocity field.


    '''
    pass


class Ellipse(namedtuple("Ellipse", ['x0', 'y0', 'a', 'b'])):
    '''
    A class that represents an ellipse.

    Attributes
    ----------

    x0 : float
    The x-coordinate of the center of the ellipse.

    y0 : float
    The y-coordinate of the center of the ellipse.

    a : float
    The major axis of the ellipse.

    b : float
    The minor axis of the ellipse.


    '''
    pass


class EllipseProperties(namedtuple("BodyProperties", ['x_cor', 'y_cor', 'u', 'v', 'circulation'])):
    '''
    A class that represents the properties of an ellipse discretized and placed in a flow field.

    Attributes
    ----------

    x_cor : np.ndarray
    The x-array coordinates of the ellipse.

    y_cor : np.ndarray
    The y-array coordinates of the ellipse.

    u : np.ndarray
    The x-component of the velocity field on the ellipse.

    v : np.ndarray
    The y-component of the velocity field on the ellipse.

    circulation : float
    The circulation of the ellipse.

    '''
    pass


@dataclass
class Geometry:
    '''
    A class that represents a geometry.

    Attributes
    ----------

    x : np.ndarray
    The x-array coordinates of the geometry.

    y : np.ndarray
    The y-array coordinates of the geometry.

    AoA : float
    The angle of attack of the geometry.
    '''
    x: np.ndarray
    y: np.ndarray
    AoA: float


class PanelizedGeometry(tp.NamedTuple('PanelizedGeometry', [('panel_length', np.ndarray),
                                                            ('panel_orientation_angle', np.ndarray),
                                                            ('panel_normal_angle', np.ndarray),
                                                            ('beta', np.ndarray),
                                                            ('control_points_x_cor', np.ndarray),
                                                            ('control_points_y_cor', np.ndarray)])):
    """
    A class that represents a panelized geometry.

    Attributes
    ----------

    panel_length : np.ndarray

    panel_orientation_angle : np.ndarray

    panel_normal_angle : np.ndarray

    beta : np.ndarray
    The angle between the freestream and the panel normal.

    control_points_x_cor : np.ndarray

    control_points_y_cor : np.ndarray

    """
    pass
