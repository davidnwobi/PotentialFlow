import pytest
import numpy as np
import random
from src import elementary_flows as ef
from src import flow_field as ff
from scipy import optimize


def find_closest_value_index(array, target):
    index = np.searchsorted(array, target)
    if index >= len(array):
        return len(array) - 1
    elif index == 0:
        return 0
    else:
        if abs(array[index] - target) < abs(array[index - 1] - target):
            return index
        else:
            return index - 1


@pytest.fixture
def default_data(request):
    num_points = 10000
    x = np.linspace(-3, 3, num=num_points)
    y = np.linspace(-3, 3, num=num_points)
    X, Y = np.meshgrid(x, y)

    yield x, y, X, Y, num_points


def test_non_lifting_flow(default_data):
    x, y, X, Y, num_points = default_data

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
    stream_function = flow.stream_function(X, Y)
    stag_x1 = find_closest_value_index(x, -radius)
    stag_x2 = find_closest_value_index(x, radius)
    stag_y = find_closest_value_index(y, 0)

    print(velocity_field[0][stag_y, stag_x1])
    assert velocity_field[0][stag_y, stag_x1] == pytest.approx(0, abs=1e-2)
    assert velocity_field[0][stag_y, stag_x2] == pytest.approx(0, abs=1e-2)
    assert stream_function[stag_y, stag_x1] == pytest.approx(0, abs=1e-2)


def test_lifting_flow(default_data):
    x, y, X, Y, num_points = default_data

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
    stream_function = flow.stream_function(X, Y)
    stag_x1 = find_closest_value_index(x, 0)
    stag_y = find_closest_value_index(y, -radius)

    assert velocity_field[0][stag_y, stag_x1] == pytest.approx(0, abs=1e-2)


def test_rankine_oval(default_data):
    x, y, X, Y, num_points = default_data

    strength = 7
    location = 1.5
    velocity = 1
    v1 = ef.Source(x_pos=-location, y_pos=0, strength=strength)
    v2 = ef.Source(x_pos=location, y_pos=0, strength=-strength)  # Negative strength is a sink
    u1 = ef.UniformFlow(horizontal_vel=velocity, vertical_vel=0)

    flow = ff.FlowField([v1, v2, u1])
    velocity_field = flow.velocity(X, Y)
    stream_function = flow.stream_function(X, Y)

    half_body_length = np.sqrt(strength / (np.pi * velocity * location) + 1)*location

    def half_body_width_func(h):
        a = location
        m = strength
        return 0.5 - (h ** 2 / a - 1) * np.tan(2 * (np.pi * velocity * a / m) * h / a) - h/a

    half_body_width = optimize.newton(half_body_width_func, half_body_length)
    loc_half_body_width = find_closest_value_index(x, half_body_length)
    zero_y = find_closest_value_index(y, 0)
    zero_x = find_closest_value_index(x, 0)
    loc_half_body_width = find_closest_value_index(y, half_body_width)

    assert velocity_field[0][zero_y, loc_half_body_width] == pytest.approx(0, abs=1e-2)
    assert stream_function[zero_y, loc_half_body_width] == pytest.approx(0, abs=1e-2)
    assert velocity_field[1][loc_half_body_width, zero_x] == pytest.approx(0, abs=1e-2)


