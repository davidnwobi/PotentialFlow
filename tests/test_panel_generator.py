import matplotlib.pyplot as plt
import copy
from src.panel_generator import PanelGenerator
from src.data_collections import PanelizedGeometry, Geometry
from src.plotting import plot_panelized_geometry
from src.airfoil_generator import generate_four_digit_NACA
import numpy as np
import pytest


@pytest.fixture()
def circle_setup_ccw():
    radius = 1
    num_panels = 9
    x = radius * np.cos(np.linspace(0, 2 * np.pi / 4, num_panels + 1))
    y = radius * np.sin(np.linspace(0, 2 * np.pi / 4, num_panels + 1))
    circle = Geometry(x, y, AoA=0)
    return circle


@pytest.fixture()
def circle_setup_cw():
    radius = 1
    num_panels = 9
    x = radius * np.flip(np.cos(np.linspace(0, 2*np.pi, num_panels + 1)))
    y = radius * np.flip(np.sin(np.linspace(0, 2*np.pi, num_panels + 1)))

    circle = Geometry(x, y, AoA=0)
    return circle


def test_is_clockwise(circle_setup_ccw):
    circle = circle_setup_ccw
    assert PanelGenerator.is_clockwise(circle) == False


def test_is_clockwise(circle_setup_cw):
    circle = circle_setup_cw
    assert PanelGenerator.is_clockwise(circle) == True


def test_ensure_cw(circle_setup_ccw):
    circle = circle_setup_ccw
    PanelGenerator.ensure_cw(circle)
    assert PanelGenerator.is_clockwise(circle) == True


def test_ensure_cw(circle_setup_cw):
    circle = circle_setup_cw
    old_circle = copy.deepcopy(circle)
    PanelGenerator.ensure_cw(circle)
    assert np.equal(circle.x, old_circle.x).all()
    assert np.equal(circle.y, old_circle.y).all()


def test_panel_generator_circle():
    panel_generator = PanelGenerator()

    # Define circle geometry
    radius = 1
    num_panels = 9
    x = radius * np.cos(np.linspace(0, 2 * np.pi, num_panels + 1)) # n+1 points for n panels
    y = radius * np.sin(np.linspace(0, 2 * np.pi, num_panels + 1))
    circle = Geometry(x, y, AoA=0)

    panelized_geometry = PanelGenerator.compute_geometric_quantities(circle)
    assert len(panelized_geometry.S) == num_panels
    assert len(panelized_geometry.phi) == num_panels
    assert len(panelized_geometry.delta) == num_panels
    assert len(panelized_geometry.beta) == num_panels
    assert len(panelized_geometry.xC) == num_panels
    assert len(panelized_geometry.yC) == num_panels

    plot_panelized_geometry(circle, panelized_geometry).show()


def test_panel_generator_airfoil():
    naca_num = 2412
    chord_length = 1
    num_panels = 100

    x, y = generate_four_digit_NACA(naca_num, num_panels+1, chord_length) # n+1 points for n panels

    airfoil = Geometry(x, y, AoA=0)

    panelized_geometry = PanelGenerator.compute_geometric_quantities(airfoil)

    assert len(panelized_geometry.S) == num_panels
    assert len(panelized_geometry.phi) == num_panels
    assert len(panelized_geometry.delta) == num_panels
    assert len(panelized_geometry.beta) == num_panels
    assert len(panelized_geometry.xC) == num_panels
    assert len(panelized_geometry.yC) == num_panels

    plot_panelized_geometry(airfoil, panelized_geometry).show()
