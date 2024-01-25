from . import geometric_integrals as gi
from . import data_collections as dc
import numpy as np
from functools import partial
from itertools import product
from multiprocessing import Pool
import numba as nb
import time


@nb.njit(cache=True)
def point_in_polygon(x: float, y: float, polygon: np.ndarray) -> bool:
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


@nb.njit(cache=True)
def compute_grid_velocity_source(panelized_geometry: dc.PanelizedGeometry, x: np.ndarray, y: np.ndarray,
                                 lam: np.ndarray,
                                 free_stream_velocity: float = 1., AoA: float = 0.) -> (np.ndarray, np.ndarray):
    Mxpj, Mypj = gi.compute_grid_geometric_integrals_source(panel_geometry=panelized_geometry, grid_x=x, grid_y=y)
    X = panelized_geometry.xC - panelized_geometry.S / 2 * np.cos(panelized_geometry.phi)
    Y = panelized_geometry.yC - panelized_geometry.S / 2 * np.sin(panelized_geometry.phi)
    np.append(X, X[0])
    np.append(Y, Y[0])
    shape = np.vstack((X, Y)).T
    u = np.zeros((len(x), len(y)))
    v = np.zeros((len(x), len(y)))
    for i in range(len(x)):
        for j in range(len(y)):
            if point_in_polygon(x=float(x[i]), y=float(y[j]), polygon=shape):
                u[j, i] = 0
                v[j, i] = 0
            else:
                u[j, i] = np.sum(lam * Mxpj[j, i] / (2 * np.pi)) + free_stream_velocity * np.cos(AoA * np.pi / 180)
                v[j, i] = np.sum(lam * Mypj[j, i] / (2 * np.pi)) + free_stream_velocity * np.sin(AoA * np.pi / 180)
    return u, v


@nb.njit(cache=True)
def compute_source_strengths(panelized_geometry: dc.PanelizedGeometry, V: float, I: np.ndarray) -> np.ndarray:
    A = I
    b = -np.cos(panelized_geometry.beta) * V * 2 * np.pi
    lam = np.linalg.solve(A, b)
    return lam


@nb.njit(cache=True)
def compute_panel_velocities_source(panelized_geometry: dc.PanelizedGeometry, lam: np.ndarray, V: float, I: np.ndarray,
                                    J: np.ndarray) -> (np.ndarray, np.ndarray):
    V_normal = np.empty(len(panelized_geometry.xC))
    V_tangential = np.empty(len(panelized_geometry.xC))
    for i in range(len(panelized_geometry.xC)):
        V_normal[i] = np.sum(lam * I[i] / (2 * np.pi)) + V * np.cos(panelized_geometry.beta[i])
        V_tangential[i] = np.sum(lam * J[i] / (2 * np.pi)) + V * np.sin(panelized_geometry.beta[i])
    return V_normal, V_tangential


def run_source_panel_method(panelized_geometry: dc.PanelizedGeometry, V: float, AoA: float, x,
                            y) -> dc.SourcePanelMethodResults:
    I, J = gi.compute_panel_geometric_integrals_source(panel_geometry=panelized_geometry)
    lam = compute_source_strengths(panelized_geometry=panelized_geometry, V=V, I=I)
    V_normal, V_tangential = compute_panel_velocities_source(panelized_geometry=panelized_geometry, lam=lam, V=V, I=I,
                                                             J=J)
    u, v = compute_grid_velocity_source(panelized_geometry=panelized_geometry, x=x, y=y, lam=lam,
                                        free_stream_velocity=V,
                                        AoA=AoA)

    panel_results = dc.SourcePanelMethodResults(V_normal=V_normal, V_tangential=V_tangential,
                                                Source_Strengths=lam, V_horizontal=u, V_vertical=v)

    return panel_results
