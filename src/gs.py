# Groud state properties, calculate c_n(0) and omega_0

import numpy as np


def delta(x, y):
    return np.where(x == y, 1, 0)


def psi0(cns, N):
    i = np.arange(1, N)
    return np.sum(np.sqrt(i) * cns[i - 1] * cns[i])


def cns(max_iter, N, dJU, muU):
    mat = np.zeros((N, N))
    psi0 = 1
    i = np.arange(N)
    j = np.arange(N)

    for k in range(max_iter):
        i_matrix, j_matrix = np.meshgrid(i, j)
        mat = (
            (0.5 * i_matrix * (i_matrix - 1) - muU * i_matrix)
            * delta(i_matrix, j_matrix)
            - dJU * np.sqrt(i_matrix) * psi0 * delta(i_matrix, j_matrix + 1)
            - dJU * np.sqrt(i_matrix + 1) * psi0 * delta(i_matrix + 1, j_matrix)
        )

        Eigvals, Eigvecs = np.linalg.eig(mat)
        ind = np.argmin(Eigvals)
        cns = np.abs(np.real(Eigvecs[:, ind]))

        psi0 = psi0(cns, N)

        # omega0U = 0
        # ns = np.arange(0, N, 1)
        # n0 =  np.sum(ns * cns * cns)
        # omega0U = np.sum( ( 0.5 * ns * (ns - 1) - muU * ns)  * cns * cns)

        # n0s[count] = n0
        # omega0U = - 2 * dJU * psi0 * psi0 + omega0U
    print("dJU = ", dJU, "psi0 = ", psi0, np.size(mat))
    return cns


def omega0U(max_iter, N, dJU, muU):
    cns = cns(max_iter, N, dJU, muU)
    psi0 = psi0(cns, N)
    ns = np.arange(0, N, 1)
    n0 = np.sum(ns * cns * cns)
    omega0U = -2 * dJU * psi0 * psi0 + np.sum(
        (0.5 * ns * (ns - 1) - muU * ns) * cns * cns
    )
    return omega0U
