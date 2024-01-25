import numpy as np
from . import data_collections as dc


class PanelGenerator:
    @staticmethod
    def is_clockwise(geometry: dc.Geometry):
        edge = (geometry.x * np.roll(geometry.y, 1)) - (geometry.y * np.roll(geometry.x, 1))  # Shoelace formula
        return np.sum(edge) > 0

    @staticmethod
    def ensure_cw(geometry: dc.Geometry):
        if not PanelGenerator.is_clockwise(geometry):
            geometry.x = np.flip(geometry.x)
            geometry.y = np.flip(geometry.y)

    @staticmethod
    def compute_geometric_quantities(geometry: dc.Geometry) -> dc.PanelizedGeometry:
        PanelGenerator.ensure_cw(geometry)
        control_points_x_cor = 0.5 * (geometry.x[:-1] + geometry.x[1:])
        control_points_y_cor = 0.5 * (geometry.y[:-1] + geometry.y[1:])
        dx = np.diff(geometry.x)
        dy = np.diff(geometry.y)

        panel_length = np.sqrt(dx ** 2 + dy ** 2)  # S
        panel_orientation_angle = np.arctan2(dy, dx)  # phi [rad] (angle between x-axis and inside panel surface)
        panel_orientation_angle = np.where(panel_orientation_angle < 0, panel_orientation_angle + 2 * np.pi,
                                           panel_orientation_angle)  # Add 2pi to the panel angle if it is negative
        panel_normal_angle = panel_orientation_angle + (
                np.pi / 2)  # delta [rad] (angle between x-axis and outside panel normal)

        beta = panel_normal_angle - (
                geometry.AoA * (np.pi / 180))  # Angle between freestream and outsidepanel normal [rad]
        beta = np.where(beta > 2 * np.pi, beta - 2 * np.pi, beta)
        beta = np.where(beta < 0, beta + 2 * np.pi, beta)

        panel_geometry = dc.PanelizedGeometry(panel_length, panel_orientation_angle, panel_normal_angle, beta,
                                              control_points_x_cor, control_points_y_cor)

        return panel_geometry
