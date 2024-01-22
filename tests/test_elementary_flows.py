import pytest
import numpy as np
import random
from src import circulation as c
from src import elementary_flows as ef
from src import data_collections as dc
from src import flow_field as ff


@pytest.fixture
def default_data(request):
    num_points = 100
    x = np.linspace(-3, 3, num=num_points)
    y = np.linspace(-3, 3, num=num_points)
    X, Y = np.meshgrid(x, y)

    yield x, y, X, Y, num_points


def test_uniform_flow_stream_fucntion(default_data):
    x, y, X, Y, num_points = default_data
    # Define the flow field
    horizontal_vel = 10
    u1 = ef.UniformFlow(horizontal_vel=horizontal_vel, vertical_vel=0)
    stream_function = np.array(u1.stream_function(X, Y))

    indices1 = np.where(stream_function == horizontal_vel)
    indices2 = np.where(y == 1)
    assert np.allclose(indices1[0], indices2)

    indices1 = np.where(stream_function == -3 * horizontal_vel)
    indices2 = np.where(y == -3)

    assert np.allclose(indices1[0], indices2)


def test_uniform_flow_velocity(default_data):
    x, y, X, Y, num_points = default_data
    # Define the flow field
    horizontal_vel = 10
    u1 = ef.UniformFlow(horizontal_vel=horizontal_vel, vertical_vel=0)
    velocity_field = np.array(u1.velocity(X, Y))  # [0] is u, [1] is v

    # Check that the horizontal velocity is constant
    assert np.allclose(velocity_field[0], horizontal_vel)

    # Check that the vertical velocity is zero
    assert np.allclose(velocity_field[1], 0)


def test_vortex_stream_fucntion(default_data):
    x, y, X, Y, num_points = default_data
    # Define the flow field
    v1 = ef.Vortex(x_pos=0, y_pos=0, circulation=10)
    stream_function = np.array(v1.stream_function(X, Y))

    # Mirror points are on the same stream function contour
    x1 = 0
    y1 = 0
    x2 = num_points - 1
    y2 = num_points - 1

    assert stream_function[x1, y1] == pytest.approx(stream_function[x2, y2], abs=1e-6)

    # Mirror points are on the same stream function contour
    x1 = random.randint(0, num_points - 1)
    y1 = random.randint(0, num_points - 1)
    x2 = num_points - 1 - x1
    y2 = num_points - 1 - y1

    assert stream_function[y1, x1] == pytest.approx(stream_function[y2, x2], abs=1e-6)


def test_vortex_velocity(default_data):
    x, y, X, Y, num_points = default_data
    # Define the flow field
    circulation = 10
    v1 = ef.Vortex(x_pos=0, y_pos=0, circulation=circulation)
    velocity_field = np.array(v1.velocity(X, Y))  # [0] is u, [1] is v

    # Random points
    x1, y1 = random.randint(0, num_points - 1), random.randint(0, num_points - 1)

    u, v = velocity_field[0, y1, x1], velocity_field[1, y1, x1]

    assert circulation / (2 * np.pi * np.sqrt(x[x1] ** 2 + y[y1] ** 2)) == pytest.approx(np.sqrt(u ** 2 + v ** 2),
                                                                                         abs=1e-6)


def test_doublet_stream_function(default_data):
    x, y, X, Y, num_points = default_data
    # Define the flow field
    kappa = 10
    d1 = ef.Doublet(x_pos=0, y_pos=0, kappa=kappa)
    stream_function = np.array(d1.stream_function(X, Y))

    # Mirror points are on the mirrored stream function contour
    x1 = 0
    y1 = 0
    x2 = num_points - 1
    y2 = num_points - 1

    assert stream_function[y1, x1] == pytest.approx(-stream_function[y2, x2], abs=1e-6)

    # Mirror points are on the mirrored stream function contour
    x1 = random.randint(0, num_points - 1)
    y1 = random.randint(0, num_points - 1)
    x2 = num_points - 1 - x1
    y2 = num_points - 1 - y1

    assert stream_function[y1, x1] == pytest.approx(-stream_function[y2, x2], abs=1e-6)


def test_doublet_velocity(default_data):
    x, y, X, Y, num_points = default_data

    # Define the flow field
    kappa = 10
    d1 = ef.Doublet(x_pos=0, y_pos=0, kappa=kappa)
    velocity_field = np.array(d1.velocity(X, Y))  # [0] is u, [1] is v

    # Mirror points have about the x axis
    x1, y1 = random.randint(0, num_points - 1), random.randint(0, num_points - 1)
    x2, y2 = x1, num_points - 1 - y1

    u1, v1 = velocity_field[0, y1, x1], velocity_field[1, y1, x1]
    u2, v2 = velocity_field[0, y2, x2], velocity_field[1, y2, x2]

    # Same x-velocity
    assert u1 == pytest.approx(u2, abs=1e-6)
    # Opposite y-velocity
    assert v1 == pytest.approx(-v2, abs=1e-6)

    # Mirror points about the x and y axis
    x1, y1 = random.randint(0, num_points - 1), random.randint(0, num_points - 1)
    x2, y2 = num_points - 1 - x1, num_points - 1 - y1

    # Same x-velocity
    u1, v1 = velocity_field[0, y1, x1], velocity_field[1, y1, x1]
    # Same y-velocity
    u2, v2 = velocity_field[0, y2, x2], velocity_field[1, y2, x2]

    assert u1 == pytest.approx(u2, abs=1e-6)
    assert v1 == pytest.approx(v2, abs=1e-6)


def test_source_stream_function(default_data):
    x, y, X, Y, num_points = default_data
    # Define the flow field
    strength = 10
    s1 = ef.Source(x_pos=0, y_pos=0, strength=strength)
    stream_function = np.array(s1.stream_function(X, Y))

    # Mirror points about the y axis have opposite stream function
    x1 = 0
    y1 = 0
    x2 = x1
    y2 = num_points - 1 - y1

    assert stream_function[y1, x1] == pytest.approx(-stream_function[y2, x2], abs=1e-6)

    # Mirror points about the y axis have opposite stream function
    x1 = random.randint(0, num_points - 1)
    y1 = random.randint(0, num_points - 1)
    x2 = x1
    y2 = num_points - 1 - y1

    assert stream_function[y1, x1] == pytest.approx(-stream_function[y2, x2], abs=1e-6)


def test_source_velocity(default_data):
    x, y, X, Y, num_points = default_data
    # Define the flow field
    strength = 10
    s1 = ef.Source(x_pos=0, y_pos=0, strength=strength)
    velocity_field = np.array(s1.velocity(X, Y))  # [0] is u, [1] is v

    x1, y1 = random.randint(0, num_points - 1), random.randint(0, num_points - 1)

    u, v = velocity_field[0, y1, x1], velocity_field[1, y1, x1]

    assert strength / (2 * np.pi * np.sqrt(x[x1] ** 2 + y[y1] ** 2)) == pytest.approx(np.sqrt(u ** 2 + v ** 2)
                                                                                      , abs=1e-6)
