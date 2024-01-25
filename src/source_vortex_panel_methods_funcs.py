from . import geometric_integrals as gi
from .source_panel_methods_funcs import point_in_polygon
import src.data_collections as dc
import numpy as np
import numba as nb


def compute_source_vortex_strengths(panelized_geometry, V, I, J, K, L):
    np.fill_diagonal(I, np.pi)
    np.fill_diagonal(J, 0)
    np.fill_diagonal(K, 0)
    np.fill_diagonal(L, 0)
    A = np.zeros((len(I[0]) + 1, len(I[0]) + 1))
    A[:-1, :-1] = I
    A[:-1, -1] = -np.sum(K, axis=1)

    kutta_condition = np.zeros(len(A[0]))
    kutta_condition[:-1] = J[0] + J[-1]
    kutta_condition[-1] = -np.sum(L[0] + L[-1]) + 2 * np.pi
    A[-1, :] = kutta_condition
    print('Shape of A:   ', A.shape)
    b = -np.cos(panelized_geometry.beta) * V * 2 * np.pi
    b = np.append(b, -V * 2 * np.pi * (np.sin(panelized_geometry.beta[0]) + np.sin(panelized_geometry.beta[-1])))
    print('Shape of b:   ', b.shape)
    solution = np.linalg.solve(A, b)
    lam = solution[:-1]
    gamma = solution[-1]
    return lam, gamma


def compute_panel_velocities_source_vortex(panelized_geometry, lam, gamma, V, I, J, K, L):
    V_normal = np.empty(len(panelized_geometry.xC))
    V_tangential = np.empty(len(panelized_geometry.xC))
    for i in range(len(panelized_geometry.xC)):
        # Note: I don't have to add the lam[i]/2 term because it is already included in the geometric integrals when I
        # filled the diagonal with pi
        V_normal[i] = -gamma * np.sum(K[i] / (2 * np.pi)) + np.sum(lam * I[i] / (2 * np.pi)) + V * np.cos(
            panelized_geometry.beta[i])
        V_tangential[i] = -gamma * np.sum(L[i] / (2 * np.pi)) + np.sum(lam * J[i] / (2 * np.pi)) + V * np.sin(
            panelized_geometry.beta[i]) + gamma / 2
    return V_normal, V_tangential



def compute_grid_velocity_source_vortex(panelized_geometry, x, y, lam, gamma, free_stream_velocity=1, AoA=0):
    Mxpj, Mypj = gi.compute_grid_geometric_integrals_source(panel_geometry=panelized_geometry, grid_x=x, grid_y=y)
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
                u[j, i] = -np.sum(gamma * Nxpj[j, i] / (2 * np.pi)) + np.sum(
                    lam * Mxpj[j, i] / (2 * np.pi)) + free_stream_velocity * np.cos(AoA * np.pi / 180)
                v[j, i] = -np.sum(gamma * Nypj[j, i] / (2 * np.pi)) + np.sum(
                    lam * Mypj[j, i] / (2 * np.pi)) + free_stream_velocity * np.sin(AoA * np.pi / 180)
    return u, v


def run_source_vortex_panel_method(panelized_geometry: dc.PanelizedGeometry, V: float, AoA: float, x,
                            y) -> dc.SourceVortexPanelMethodResults:
    I, J = gi.compute_panel_geometric_integrals_source(panelized_geometry)
    K, L = gi.compute_panel_geometric_integrals_vortex(panelized_geometry)
    lam, gamma = compute_source_vortex_strengths(panelized_geometry, V, I, J, K, L)
    V_normal, V_tangential = compute_panel_velocities_source_vortex(panelized_geometry, lam, gamma, V, I, J, K, L)

    u, v = compute_grid_velocity_source_vortex(panelized_geometry, x, y, lam, gamma, V, AoA)

    panel_results = dc.SourceVortexPanelMethodResults(V_normal=V_normal, V_tangential=V_tangential,
                                                      Source_Strengths=lam, Circulation=gamma, V_horizontal=u,
                                                      V_vertical=v)

    return panel_results
