from . import geometric_integrals as gi
import numpy as np


def point_in_polygon(x, y, polygon):
    n = len(polygon)
    inside = False

    x1, y1 = polygon[0]
    for i in range(n + 1):
        x2, y2 = polygon[i % n]
        if y > min(y1, y2):
            if y <= max(y1, y2):
                if x <= max(x1, x2):
                    if y1 != y2:
                        x_intersect = (y - y1) * (x2 - x1) / (y2 - y1) + x1
                        if x1 == x2 or x <= x_intersect:
                            inside = not inside
        x1, y1 = x2, y2

    return inside


def compute_grid_velocity(panelized_geometry, x, y, gamma, free_stream_velocity=1, AoA=0):
    Mxpj, Mypj = gi.compute_grid_geometric_integrals(panelized_geometry, x, y)
    X = panelized_geometry.xC - panelized_geometry.S / 2 * np.cos(panelized_geometry.phi)
    Y = panelized_geometry.yC - panelized_geometry.S / 2 * np.sin(panelized_geometry.phi)
    np.append(X, X[0])
    np.append(Y, Y[0])
    shape = np.array([X, Y]).T
    u = np.zeros((len(x), len(y)))
    v = np.zeros((len(x), len(y)))
    for i in range(len(x)):
        for j in range(len(y)):
            if point_in_polygon(x[i], y[j], shape):
                u[j, i] = 0
                v[j, i] = 0
                continue
            u[j, i] = np.sum(gamma * Mxpj[j, i] / (2 * np.pi)) + free_stream_velocity * np.cos(AoA * np.pi / 180)
            v[j, i] = np.sum(gamma * Mypj[j, i] / (2 * np.pi)) + free_stream_velocity * np.sin(AoA * np.pi / 180)
    return u, v


def compute_source_strengths(panelized_geometry, V, I):
    A = I
    b = -np.cos(panelized_geometry.beta) * V * 2 * np.pi
    gamma = np.linalg.solve(A, b)
    return gamma


def compute_panel_velocities(panelized_geometry, gamma, V, I, J):
    V_normal = np.empty(len(panelized_geometry.xC))
    for i in range(len(panelized_geometry.xC)):
        V_normal[i] = np.sum(gamma * I[i] / (2 * np.pi)) + V * np.cos(panelized_geometry.beta[i])
    V_tangential = np.empty(len(panelized_geometry.xC))
    for i in range(len(panelized_geometry.xC)):
        V_tangential[i] = np.sum(gamma * J[i] / (2 * np.pi)) + V * np.sin(panelized_geometry.beta[i])
    return V_normal, V_tangential
