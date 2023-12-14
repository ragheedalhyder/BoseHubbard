# A class that includes all the groundstate functions

import numpy as np


class excitations:
    def __init__(self, grid, groundstate, cns):
        self.grid = grid
        self.groundstate = groundstate
        self.cns = cns

    # def __str__(self):
    #     return f"grid = {self.grid}, groundstate = {self.groundstate}"

    # def __repr__(self):
    #     return f"grid = {self.grid}, groundstate = {self.groundstate}"

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
                # print( x, B[5][5])

                AB = np.concatenate((A, B), axis=1)
                BA = np.concatenate((-B, -A), axis=1)
                MatAB = np.concatenate((AB, BA), axis=0)
                Eigvals, Eigvecs = np.linalg.eig(MatAB)

                Eigsorted = np.sort(Eigvals)
                omega0 = Eigsorted[N : 2 * N]
                # print(omega0)
                for lambda_ in range(N - 1):
                    omegaklambda[lambda_][kx][ky] = np.real(omega0[lambda_])
                    # print(omegaklambda[lambda_][kx][ky])
                    doublezero = False
                    if lambda_ == 0:
                        ind = np.where(
                            np.floor(10 * np.abs(np.real(Eigvals)))
                            == np.floor(10 * np.abs(omega0[lambda_]))
                        )

                        if ind[0].size > 2:
                            doublezero = True

                        for n in range(N):
                            uks[n][kx][ky][lambda_] = cns[n]
                            vks[n][kx][ky][lambda_] = 0
                        omegaklambda[lambda_][kx][ky] = 0
                    else:
                        ind = np.where(Eigvals == omega0[lambda_])
                        # print(omega0[lambda_])
                        ind1 = ind[0][0]
                        # print(ind1)
                        if ind[0].size > 1:
                            ind1 = ind[0][1]
                        Norm = 0
                        Norm = sum(
                            np.real(Eigvecs[0 : N - 1, ind1])
                            * np.real(Eigvecs[0 : N - 1, ind1])
                        ) - sum(
                            np.real(Eigvecs[N : 2 * N - 1, ind1])
                            * np.real(Eigvecs[N : 2 * N - 1, ind1])
                        )
                        # if(lambda_ == 8):
                        #     print(Eigvecs[0 : 2 * N - 1, ind1])
                        if round(Norm, 6) <= 0:
                            # print(dJU, kx, ky, lambda_, Norm, omegaklambda[lambda_][kx][ky])
                            # print(Eigvecs[0 : 2 * N - 1, ind1])
                            Norm = 1
                        uks[:][kx][ky][lambda_] = np.real(Eigvecs[0:N, ind1]) / np.sqrt(
                            Norm
                        )
                        vks[:][kx][ky][lambda_] = np.real(
                            Eigvecs[N : 2 * N, ind1]
                        ) / np.sqrt(Norm)
        return uks, vks, omegaklambda


# plt.plot(dJUs, omega0s)
# plt.plot(dJUs, omega1s)
# plt.plot(dJUs, omega2s)
# plt.show()
