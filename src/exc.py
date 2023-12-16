# A class that includes all the groundstate functions

import numpy as np


class excitations:
    def __init__(self, grid, groundstate, cns):
        self.grid = grid
        self.groundstate = groundstate
        self.cns = cns

    def __str__(self):
        return f"grid = {self.grid}, groundstate = {self.groundstate}"

    def __repr__(self):
        return f"grid = {self.grid}, groundstate = {self.groundstate}"

    def cn(self, n):
        N = self.groundstate.N
        cns = self.cns
        if n == -1:
            return 0
        elif n == N:
            return 0
        else:
            return cns[n].real

    def JkU(self, x):
        dJU = self.groundstate.dJU
        return dJU - dJU * x

    def epsI(self, kx, ky):
        return pow(np.sin(kx / 2), 2) + pow(np.sin(ky / 2), 2)

    def calculate_matrices(self):
        Lx = self.grid.Lx
        Ly = self.grid.Ly
        KXs = self.grid.KXs
        KYs = self.grid.KYs
        N = self.groundstate.N
        muU = self.groundstate.muU
        cns = self.cns

        psi0 = self.groundstate.psi0(cns)
        omega0U = self.groundstate.omega0U(cns)

        L = Lx = Ly
        A = np.zeros((N, N))
        B = np.zeros((N, N))
        uks = np.zeros((N, L, L, N))
        vks = np.zeros((N, L, L, N))
        omegaklambda = np.zeros((N, L, L))

        for kx in range(Lx):
            for ky in range(Ly):
                x = self.epsI(KXs[kx], KYs[ky])
                for n in range(N):
                    for m in range(N):
                        A[n, m] = (
                            (0.5 * n * (n - 1) - muU * n - omega0U) * (n == m)
                            - self.JkU(0)
                            * psi0
                            * (np.sqrt(n) * (n == m + 1) + np.sqrt(m) * (n + 1 == m))
                            - self.JkU(x)
                            * (
                                np.sqrt(n)
                                * np.sqrt(m)
                                * self.cn(m - 1)
                                * self.cn(n - 1)
                                + np.sqrt(n + 1)
                                * np.sqrt(m + 1)
                                * self.cn(m + 1)
                                * self.cn(n + 1)
                            )
                        )
                        B[n, m] = -self.JkU(x) * (
                            np.sqrt(n)
                            * np.sqrt(m + 1)
                            * self.cn(m + 1)
                            * self.cn(n - 1)
                            + np.sqrt(n + 1)
                            * np.sqrt(m)
                            * self.cn(m - 1)
                            * self.cn(n + 1)
                        )

                AB = np.concatenate((A, B), axis=1)
                BA = np.concatenate((-B, -A), axis=1)
                MatAB = np.concatenate((AB, BA), axis=0)
                Eigvals, Eigvecs = np.linalg.eig(MatAB)

                Eigsorted = np.sort(Eigvals)
                sorted_indices = np.argsort(Eigvals)
                omega0 = Eigsorted[N : 2 * N]
                ind = sorted_indices[N : 2 * N]

                uks[:][kx][ky][0] = cns
                vks[:][kx][ky][0] = np.zeros(N)
                omegaklambda[0][kx][ky] = 0

                for lambda_ in range(1, N - 1):
                    ind1 = ind[lambda_]
                    omegaklambda[lambda_][kx][ky] = np.real(omega0[lambda_])
                    uks_iter = np.real(Eigvecs[0:N, ind1])
                    vks_iter = np.real(Eigvecs[N : 2 * N, ind1])

                    Norm = np.dot(uks_iter, uks_iter) - np.dot(vks_iter, vks_iter)

                    if round(Norm, 6) <= 0:
                        print("Norm is negative")
                        Norm = 1

                    uks[:][kx][ky][lambda_] = uks_iter / np.sqrt(Norm)
                    vks[:][kx][ky][lambda_] = vks_iter / np.sqrt(Norm)

        return uks, vks, omegaklambda
