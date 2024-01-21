from collections import namedtuple


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
