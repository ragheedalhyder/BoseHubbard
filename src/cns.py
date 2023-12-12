import numpy as np


def delta(x, y):
    return np.where(x == y, 1, 0)


def Psi0(cns, N):
    i = np.arange(1, N)
    return np.sum(np.sqrt(i) * cns[i - 1] * cns[i])


def psi(max_iter, N, dJU, muU):
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

        psi0 = Psi0(cns, N)
    return cns
