# A class that includes all the groundstate functions

import numpy as np

class groundstate:
    def __init__(self, params):
        self.N = params.N
        self.dJU = params.dJU
        self.muU = params.muU

    def __str__(self):
        return f"N = {self.N}, dJU = {self.dJU}, muU = {self.muU}, psi0 = {self.psi0}, n0 = {self.n0}, omega0U = {self.omega0U}"

    def __repr__(self):
        return f" N = {self.N}, dJU = {self.dJU}, muU = {self.muU}, psi0 = {self.psi0}, n0 = {self.n0}, omega0U = {self.omega0U}"

    def delta(self, x, y):
        return np.where(x == y, 1, 0)

    def psi0(self, cns):
        N = self.N
        i = np.arange(1, N)
        psi0 = np.sum(np.sqrt(i) * cns[i - 1] * cns[i])
        return psi0

    def n0(self, cns):
        N = self.N
        i = np.arange(0, N)
        return np.sum(i * cns * cns)

    def omega0U(self, cns):
        N = self.N
        ns = np.arange(0, N, 1)
        psi0 = self.psi0(cns)
        return -2 * self.dJU * psi0 * psi0 + np.sum(
            (0.5 * ns * (ns - 1) - self.muU * ns) * cns * cns
        )

    def cns(self):
        N = self.N
        dJU = self.dJU
        muU = self.muU

        mat = np.zeros((N, N))
        psi0 = 1
        i = np.arange(N)
        j = np.arange(N)

        for k in range(20000):
            i_matrix, j_matrix = np.meshgrid(i, j)
            mat = (
                (0.5 * i_matrix * (i_matrix - 1) - muU * i_matrix)
                * self.delta(i_matrix, j_matrix)
                - dJU * np.sqrt(i_matrix) * psi0 * self.delta(i_matrix, j_matrix + 1)
                - dJU
                * np.sqrt(i_matrix + 1)
                * psi0
                * self.delta(i_matrix + 1, j_matrix)
            )

            Eigvals, Eigvecs = np.linalg.eig(mat)
            ind = np.argmin(Eigvals)
            cns = np.abs(np.real(Eigvecs[:, ind]))
            psi0 = self.psi0(cns)
        return cns
