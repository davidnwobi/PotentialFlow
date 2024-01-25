from collections import namedtuple
import numpy as np
from dataclasses import dataclass
import typing as tp
from numba import float64, int64, boolean
from numba.experimental import jitclass

SourcePanelMethodResults = namedtuple('PanelMethodResults',
                                      ['V_normal', 'V_tangential', 'Source_Strengths', 'V_horizontal', 'V_vertical'])

VortexPanelMethodResults = namedtuple('PanelMethodResults',
                                      ['V_normal', 'V_tangential', 'Vortex_Strengths', 'V_horizontal', 'V_vertical'])


class FlowFieldProperties(namedtuple("FlowField", ['x', 'y', 'u', 'v'])):
    """
    A class that represents a flow field properties.
    For x and y, DO NOT PASS IN MESHGRIDDED ARRAYS. Pass in the x and y arrays themselves.
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


    """
    pass


class Ellipse(namedtuple("Ellipse", ['x0', 'y0', 'a', 'b'])):
    """
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


    """
    pass


class EllipseProperties(namedtuple("BodyProperties", ['x_cor', 'y_cor', 'u', 'v', 'circulation'])):
    """
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

    """
    pass


specGeometry = [
    ('x', float64[:]),
    ('y', float64[:]),
    ('AoA', float64),
]
jitclass(specGeometry)


class Geometry:
    """
    A class that represents a geometry.

    Attributes
    ----------

    x : np.ndarray
    The x-array coordinates of the geometry.

    y : np.ndarray
    The y-array coordinates of the geometry.

    AoA : float
    The angle of attack of the geometry.
    """
    x: np.ndarray
    y: np.ndarray
    AoA: float

    def __init__(self, x, y, AoA):
        self.x = x
        self.y = y
        self.AoA = AoA


specPanelizedGeometry = [
    ('S', float64[:]),
    ('phi', float64[:]),
    ('delta', float64[:]),
    ('beta', float64[:]),
    ('xC', float64[:]),
    ('yC', float64[:]),
]


@jitclass(specPanelizedGeometry)
class PanelizedGeometry():
    """
    A class that represents a panelized geometry.

    Attributes
    ----------

    S : np.ndarray
    The length of the panel.

    phi : np.ndarray
    The angle between the x-axis and the inside panel surface.

    delta : np.ndarray
    The angle between the x-axis and the outside panel normal.

    beta : np.ndarray
    The angle between the freestream and the panel normal.

    xC : np.ndarray
    x-coordinate of the control point.

    yC : np.ndarray
    y-coordinate of the control point.

    """
    S: np.ndarray
    phi: np.ndarray
    delta: np.ndarray
    beta: np.ndarray
    xC: np.ndarray
    yC: np.ndarray

    def __init__(self, S, phi, delta, beta, xC, yC):
        self.S = S
        self.phi = phi
        self.delta = delta
        self.beta = beta
        self.xC = xC
        self.yC = yC
