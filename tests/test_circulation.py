import numpy as np
import pytest
from src import circulation as c
from src import elementary_flows as ef
from src import data_collections as dc
from src import flow_field as ff


@pytest.fixture
def default_data(request):
    x = np.linspace(-3, 3, num=100)
    y = np.linspace(-3, 3, num=100)
    X, Y = np.meshgrid(x, y)

    yield x, y, X, Y


def test_circulation_uniform_flow(default_data):
    """
    Tests the circulation of an ellipse placed in a uniform flow field.
    """

    x, y, X, Y = default_data
    # Define the flow field
    u1 = ef.UniformFlow(horizontal_vel=1, vertical_vel=0)
    velocity_field = u1.velocity(X, Y)

    # Define the ellipse
    ellipse_def = dc.Ellipse(x0=0, y0=0, a=1, b=1)

    flow_properties = dc.FlowFieldProperties(x, y, velocity_field[0], velocity_field[1])
    # Compute the ellipse properties
    ellipse_properties = c.compute_ellipse_and_circulation(flow_properties, ellipse_def)

    # Assert the circulation is zero
    assert ellipse_properties.circulation == pytest.approx(0, abs=1e-6)


def test_circulation_non_lifting_flow(default_data):
    """
    Tests the circulation of an ellipse placed in a non-lifting flow field.
    """
    x, y, X, Y = default_data

    velocity = 10
    radius = 1
    kappa = 2 * np.pi * velocity * radius ** 2  # Known solution for a cylinder
    vortex_strength = 4 * np.pi * velocity * radius  # Known solution for a cylinder

    alpha = 0 * np.pi / 30  # Angle of attack
    u1 = ef.UniformFlow(horizontal_vel=velocity * np.cos(alpha), vertical_vel=velocity * np.sin(alpha))
    v1 = ef.Vortex(x_pos=0, y_pos=0, circulation=vortex_strength)
    d1 = ef.Doublet(x_pos=0, y_pos=0, kappa=kappa)

    flow = ff.FlowField([d1, u1])
    velocity_field = flow.velocity(X, Y)

    # Create ellipse and compute circulation
    flow_properties = dc.FlowFieldProperties(x, y, velocity_field[0], velocity_field[1])
    ellipse = dc.Ellipse(x0=0, y0=0, a=1, b=1)
    ellipse_properties = c.compute_ellipse_and_circulation(flow_properties, ellipse, divsions=1000)

    # Assert the circulation is zero

    assert ellipse_properties.circulation == pytest.approx(0, abs=1e-6)


def test_circulation_lifting_flow(default_data):
    """
    Tests the circulation of an ellipse placed in a lifting flow field.
    """
    x, y, X, Y = default_data
    #
    velocity = 10
    radius = 1
    kappa = 2 * np.pi * velocity * radius ** 2  # Known solution for a cylinder
    vortex_strength = 4 * np.pi * velocity * radius  # Known solution for a cylinder

    alpha = 0 * np.pi / 30  # Angle of attack
    u1 = ef.UniformFlow(horizontal_vel=velocity * np.cos(alpha), vertical_vel=velocity * np.sin(alpha))
    v1 = ef.Vortex(x_pos=0, y_pos=0, circulation=vortex_strength)
    d1 = ef.Doublet(x_pos=0, y_pos=0, kappa=kappa)

    flow = ff.FlowField([d1, u1, v1])
    velocity_field = flow.velocity(X, Y)

    # Create ellipse and compute circulation
    flow_properties = dc.FlowFieldProperties(x, y, velocity_field[0], velocity_field[1])
    ellipse = dc.Ellipse(x0=0, y0=0, a=1, b=1)
    ellipse_properties = c.compute_ellipse_and_circulation(flow_properties, ellipse, divsions=1000)

    # Circulation is equal to the vortex strength when enclosing the vortex

    assert ellipse_properties.circulation == pytest.approx(vortex_strength, abs=1e-2)

    # Circulation equal to zero when not enclosing the vortex

    ellipse = dc.Ellipse(x0=-2, y0=-2, a=1, b=1)
    ellipse_properties = c.compute_ellipse_and_circulation(flow_properties, ellipse, divsions=1000)

    assert ellipse_properties.circulation == pytest.approx(0, abs=1e-6)
