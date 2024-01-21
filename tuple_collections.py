from collections import namedtuple

FlowFieldProperties = namedtuple("FlowField", ['x', 'y', 'u', 'v'])
Ellipse = namedtuple("Ellipse", ['x0', 'y0', 'a', 'b'])
EllipseProperties = namedtuple("BodyProperties", ['x_cor', 'y_cor', 'u', 'v', 'circulation'])