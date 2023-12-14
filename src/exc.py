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

    def calculate_matrices(self, muU, omega0U, dJU, psi0, A, B, cns, epsI, JkU, cn):
        Lx = self.grid.Lx
        Ly = self.grid.Ly
        KXs = self.grid.KXs
        KYs = self.grid.KYs
        muU = self.groundstate.muU
        cns = self.groundstate.cns

        N = self.groundstate.N
        muU = self.groundstate.muU

        L = Lx = Ly
        A = np.zeros((N, N))
        B = np.zeros((N, N))
        uks = np.zeros((N, L, L, N))
        vks = np.zeros((N, L, L, N))
        omegaklambda = np.zeros((N, L, L))

        for kx in range(Lx):
            for ky in range(Ly):
                x = epsI(KXs[kx], KYs[ky])
                for n in range(N):
                    for m in range(N):
                        A[n, m] = (
                            (0.5 * n * (n - 1) - muU * n - omega0U) * (n == m)
                            - JkU(dJU, 0)
                            * psi0
                            * (np.sqrt(n) * (n == m + 1) + np.sqrt(m) * (n + 1 == m))
                            - JkU(dJU, x)
                            * (
                                np.sqrt(n) * np.sqrt(m) * cn(m - 1) * cn(n - 1)
                                + np.sqrt(n + 1)
                                * np.sqrt(m + 1)
                                * cn(m + 1)
                                * cn(n + 1)
                            )
                        )
                        B[n, m] = -JkU(dJU, x) * (
                            np.sqrt(n) * np.sqrt(m + 1) * cn(m + 1) * cn(n - 1)
                            + np.sqrt(n + 1) * np.sqrt(m) * cn(m - 1) * cn(n + 1)
                        )

                AB = np.concatenate((A, B), axis=1)
                BA = np.concatenate((-B, -A), axis=1)
                MatAB = np.concatenate((AB, BA), axis=0)
                Eigvals, Eigvecs = np.linalg.eig(MatAB)

                Eigsorted = np.sort(Eigvals)
                omega0 = Eigsorted[N : 2 * N]

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

        omega0s[count] = omegaklambda[1][5][5]
        omega1s[count] = omegaklambda[2][5][5]
        omega2s[count] = omegaklambda[3][5][5]
        return uks, vks, omegaklambda


dkx = dky = 2 * np.pi / L

KXs = KYs = np.arange(-np.pi, np.pi, dkx)
dkxs = dkys = np.ones(L) * dkx
# print(KXs)

# Trapezoidal rule
dkxs[0] = dkxs[L - 1] = dkx / 2
dkys[0] = dkys[L - 1] = dky / 2

dJUs = np.arange(0.01, 0.6, 0.01)
muU = np.sqrt(2) - 1


n0s = np.zeros(len(dJUs))
omega0s = np.zeros(len(dJUs))
omega1s = np.zeros(len(dJUs))
omega2s = np.zeros(len(dJUs))

for count in range(len(dJUs)):  # (0, len(dJUs)):
    dJU = dJUs[count]
    mat = np.zeros((N, N))
    A = np.zeros((N, N))
    B = np.zeros((N, N))
    uks = np.zeros((N, L, L, N))
    vks = np.zeros((N, L, L, N))
    omegaklambda = np.zeros((N, L, L))

    # psi0 = 1

    # for k in range(5000):
    #     for i in range(N):
    #         for j in range(N):
    #             mat[i,j] = (0.5 * i * (i - 1) - muU * i) * delta(i,j) -  dJU * np.sqrt(i) * psi0 * delta(i,j + 1) - dJU * np.sqrt(i + 1) * psi0 * (i + 1 == j)

    #     Eigvals, Eigvecs = np.linalg.eig(mat)
    #     ind = np.argmin(Eigvals)
    #     cns = np.abs(np.real(Eigvecs[:,ind]))

    #     psi0 = Psi0(cns)
    # omega0U = 0
    # ns = np.arange(0, N, 1)
    # n0 =  np.sum(ns * cns * cns)
    # omega0U = np.sum( ( 0.5 * ns * (ns - 1) - muU * ns)  * cns * cns)

    # n0s[count] = n0
    # omega0U = - 2 * dJU * psi0 * psi0 + omega0U

    for kx in range(L):
        for ky in range(L):
            x = epsI(KXs[kx], KYs[ky])
            for n in range(N):
                for m in range(N):
                    A[n, m] = (
                        (0.5 * n * (n - 1) - muU * n - omega0U) * (n == m)
                        - JkU(dJU, 0)
                        * psi0
                        * (np.sqrt(n) * (n == m + 1) + np.sqrt(m) * (n + 1 == m))
                        - JkU(dJU, x)
                        * (
                            np.sqrt(n) * np.sqrt(m) * cn(m - 1) * cn(n - 1)
                            + np.sqrt(n + 1) * np.sqrt(m + 1) * cn(m + 1) * cn(n + 1)
                        )
                    )
                    B[n, m] = -JkU(dJU, x) * (
                        np.sqrt(n) * np.sqrt(m + 1) * cn(m + 1) * cn(n - 1)
                        + np.sqrt(n + 1) * np.sqrt(m) * cn(m - 1) * cn(n + 1)
                    )

            AB = np.concatenate((A, B), axis=1)
            BA = np.concatenate((-B, -A), axis=1)
            MatAB = np.concatenate((AB, BA), axis=0)
            Eigvals, Eigvecs = np.linalg.eig(MatAB)
            # sort the positive eigenvalues

            omega00 = np.sort(Eigvals)
            omega0 = omega00[N : 2 * N]

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

        omega0s[count] = omegaklambda[1][5][5]
        omega1s[count] = omegaklambda[2][5][5]
        omega2s[count] = omegaklambda[3][5][5]

plt.plot(dJUs, omega0s)
plt.plot(dJUs, omega1s)
plt.plot(dJUs, omega2s)
plt.show()
