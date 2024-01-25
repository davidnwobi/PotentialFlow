from . import geometric_integrals as gi
from .source_panel_methods_funcs import point_in_polygon
import src.data_collections as dc
import numpy as np
import numba as nb


@nb.njit(cache=True)
def compute_vortex_strengths(panelized_geometry, V, K):
    A = -K
    kutta_condition = np.zeros(len(K[0]))
    kutta_condition[0] = 1
    kutta_condition[-1] = 1
    A[-1, :] = kutta_condition
    print('A:   ', A)
    b = -np.cos(panelized_geometry.beta) * V * 2 * np.pi
    b[-1] = 0
    gamma = np.linalg.solve(A, b)
    return gamma


@nb.njit(cache=True)
def compute_grid_velocity_vortex(panelized_geometry, x, y, gamma, free_stream_velocity=1, AoA=0):
    Nxpj, Nypj = gi.compute_grid_geometric_integrals_vortex(panelized_geometry, x, y)
    X = panelized_geometry.xC - panelized_geometry.S / 2 * np.cos(panelized_geometry.phi)
    Y = panelized_geometry.yC - panelized_geometry.S / 2 * np.sin(panelized_geometry.phi)
    X = np.append(X, X[0])
    Y = np.append(Y, Y[0])
    shape = np.vstack((X, Y)).T
    u = np.zeros((len(x), len(y)))
    v = np.zeros((len(x), len(y)))
    for i in range(len(x)):
        for j in range(len(y)):
            if point_in_polygon(float(x[i]), float(y[j]), shape):
                u[j, i] = 0
                v[j, i] = 0
            else:
                u[j, i] = -np.sum(gamma * Nxpj[j, i] / (2 * np.pi)) + free_stream_velocity * np.cos(AoA * np.pi / 180)
                v[j, i] = -np.sum(gamma * Nypj[j, i] / (2 * np.pi)) + free_stream_velocity * np.sin(AoA * np.pi / 180)
    return u, v


@nb.njit(cache=True)
def compute_panel_velocities_vortex(panelized_geometry, gamma, V, K, L):
    V_normal = np.empty(len(panelized_geometry.xC))
    V_tangential = np.empty(len(panelized_geometry.xC))
    for i in range(len(panelized_geometry.xC)):
        V_normal[i] = np.sum(-gamma * K[i] / (2 * np.pi)) + V * np.cos(panelized_geometry.beta[i])
        V_tangential[i] = np.sum(-gamma * L[i] / (2 * np.pi)) + V * np.sin(panelized_geometry.beta[i]) + gamma[i] / 2

    return V_normal, V_tangential


def run_vortex_panel_method(panelized_geometry: dc.PanelizedGeometry, V: float, AoA: float, x,
                            y) -> dc.VortexPanelMethodResults:
    K, L = gi.compute_panel_geometric_integrals_vortex(panelized_geometry)
    gamma = compute_vortex_strengths(panelized_geometry, V, K)
    V_normal, V_tangential = compute_panel_velocities_vortex(panelized_geometry, gamma, V, K, L)

    u, v = compute_grid_velocity_vortex(panelized_geometry, x, y, gamma, V, AoA)

    panel_results = dc.VortexPanelMethodResults(V_normal=V_normal, V_tangential=V_tangential,
                                                Vortex_Strengths=gamma, V_horizontal=u, V_vertical=v)

    return panel_results
