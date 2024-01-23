import numpy as np
from . import data_collections as dc


def compute_geometric_integral(A, B, C, D, E, S):
    I = C * 0.5 * (np.log((S ** 2 + 2 * A * S + B) / B)) + (D - A * C) / E * (
            np.arctan2((S + A), E) - np.arctan2(A, E))
    return I


def compute_repeating_terms(x_i, y_i, X_j, Y_j, phi_j):
    A = -(x_i - X_j) * np.cos(phi_j) - (y_i - Y_j) * np.sin(phi_j)
    B = (x_i - X_j) ** 2 + (y_i - Y_j) ** 2
    E = np.sqrt(B - A ** 2)
    E = np.where(np.iscomplex(E) | np.isnan(E) | np.isinf(E), 0, E)
    return A, B, E


def normal_geometric_integral(x_i, y_i, X_j, Y_j, phi_i, phi_j, S_j, A, B, E):
    Cn = np.sin(phi_i - phi_j)
    Dn = -(x_i - X_j) * np.sin(phi_j) + (y_i - Y_j) * np.cos(phi_j)
    I_ij = compute_geometric_integral(A, B, Cn, Dn, E, S_j)

    return I_ij


def tangential_geometric_integral(x_i, y_i, X_j, Y_j, phi_i, phi_j, S_j, A, B, E):
    Ct = -np.cos(phi_i - phi_j)
    Dt = (x_i - X_j) * np.cos(phi_j) + (y_i - Y_j) * np.sin(phi_j)
    J_ij = compute_geometric_integral(A, B, Ct, Dt, E, S_j)

    return J_ij


def horizontal_geometric_integral(x_p, X_j, phi_j, S_j, A, B, E):
    Cx = -np.cos(phi_j)
    Dx = (x_p - X_j)
    M_xpj = compute_geometric_integral(A, B, Cx, Dx, E, S_j)
    return M_xpj


def vertical_geometric_integral(y_p, Y_j, phi_j, S_j, A, B, E):
    Cy = -np.sin(phi_j)
    Dy = (y_p - Y_j)
    M_ypj = compute_geometric_integral(A, B, Cy, Dy, E, S_j)
    return M_ypj


def compute_panel_geometric_integrals(panel_geometry: dc.PanelizedGeometry):
    I = np.empty((panel_geometry.S.size, panel_geometry.S.size))
    J = np.empty((panel_geometry.S.size, panel_geometry.S.size))
    X = panel_geometry.xC - panel_geometry.S / 2 * np.cos(panel_geometry.phi)
    Y = panel_geometry.yC - panel_geometry.S / 2 * np.sin(panel_geometry.phi)
    for i in range(panel_geometry.S.size):
        A, B, E = compute_repeating_terms(panel_geometry.xC[i], panel_geometry.yC[i], X,
                                          Y, panel_geometry.phi)
        null_index = np.where(np.isnan(E) | np.isinf(E) | np.iscomplex(E))
        I[i] = normal_geometric_integral(panel_geometry.xC[i], panel_geometry.yC[i], X,
                                         Y, panel_geometry.phi[i], panel_geometry.phi, panel_geometry.S,
                                         A, B, E)

        J[i] = tangential_geometric_integral(panel_geometry.xC[i], panel_geometry.yC[i], X,
                                             Y, panel_geometry.phi[i], panel_geometry.phi,
                                             panel_geometry.S,
                                             A, B, E)

        I[i, null_index] = 0
        J[i, null_index] = 0

    I = np.where(np.isnan(I) | np.isinf(I) | np.iscomplex(I), 0, I)
    J = np.where(np.isnan(J) | np.isinf(J) | np.iscomplex(J), 0, J)
    return I, J


def compute_grid_geometric_integrals(panel_geometry: dc.PanelizedGeometry, grid_x: np.ndarray, grid_y: np.ndarray):
    Ixpj = np.empty((grid_y.size, grid_x.size, panel_geometry.S.size))
    Jypj = np.empty((grid_y.size, grid_x.size, panel_geometry.S.size))
    X = panel_geometry.xC - panel_geometry.S / 2 * np.cos(panel_geometry.phi)
    Y = panel_geometry.yC - panel_geometry.S / 2 * np.sin(panel_geometry.phi)
    for i in range(grid_x.size):
        for j in range(grid_y.size):
            A, B, E = compute_repeating_terms(grid_x[i], grid_y[j], X,
                                              Y, panel_geometry.phi)

            M_xpj = horizontal_geometric_integral(grid_x[i], X, panel_geometry.phi, panel_geometry.S,
                                                  A, B, E)
            M_ypj = vertical_geometric_integral(grid_y[j], Y, panel_geometry.phi, panel_geometry.S,
                                                A, B, E)

            M_xpj[np.isnan(M_xpj)] = 0
            M_ypj[np.isnan(M_ypj)] = 0
            Ixpj[j, i], Jypj[j, i] = M_xpj, M_ypj

    Ixpj = np.where(np.isnan(Ixpj) | np.isinf(Ixpj) | np.iscomplex(Ixpj), 0, Ixpj)
    Jypj = np.where(np.isnan(Jypj) | np.isinf(Jypj) | np.iscomplex(Jypj), 0, Jypj)
    return Ixpj, Jypj
