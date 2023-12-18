# Perturbative diagrams
import numpy as np


class perturbative:
    def __init__(self, grid, params, vertices, omegaklambda):
        self.grid = grid
        self.Lx = grid.Lx
        self.Ly = grid.Ly
        self.N = params.N
        self.dJU = params.dJU
        self.muU = params.muU
        self.UIB = params.UIB
        self.cutoff = params.cutoff
        self.vertices = vertices
        self.omegaklambda = omegaklambda

    def __str__(self):
        return f"groundstate = {self.groundstate}, UIB = {self.UIB}, cutoff = {self.cutoff}"

    def __repr__(self):
        return f"groundstate = {self.groundstate}, UIB = {self.UIB}, cutoff = {self.cutoff}"

    def epsI(self, kx, ky):
        return pow(np.sin(kx / 2), 2) + pow(np.sin(ky / 2), 2)

    # Sigma0 zeroth term of self energy
    def sigma0(self, n0):
        deltan2 = 0
        M = self.grid.M
        cutoff = self.cutoff
        Lx = self.Lx
        Ly = self.Ly
        UIB = self.UIB
        cutoff = self.cutoff

        for lambda_ in range(1, cutoff):
            for kx in range(Lx):
                # print("kx = ", KXs[kx])
                for ky in range(Ly):
                    V_vert = self.vertices.V(kx, ky, lambda_, kx, ky, lambda_)
                    deltan2 = deltan2 + V_vert
        return UIB * (n0 + deltan2 / M)

    def sigma1(self):
        cutoff = self.cutoff
        Lx = self.Lx
        Ly = self.Ly
        UIB = self.UIB
        M = self.grid.M
        omegaklambda = self.omegaklambda
        cutoff = self.cutoff
        KXs = self.grid.KXs
        KYs = self.grid.KYs
        dJU = self.dJU
        Sigma1 = 0

        for lambda_ in range(1, cutoff):
            for kx in range(Lx):
                for ky in range(Ly):
                    Nkres = self.vertices.Nk(kx, ky, lambda_)
                    eps = self.epsI(KXs[kx], KYs[ky])
                    res = np.real(
                        pow(UIB, 2)
                        * abs(pow(Nkres, 2))
                        / (-omegaklambda[lambda_][kx][ky] - dJU * eps + 1j * 0.0001)
                    )
                    Sigma1 = Sigma1 + res
        return Sigma1 / M

    def sigma2(self):
        cutoff = self.cutoff
        Lx = self.Lx
        Ly = self.Ly
        M = self.grid.M
        UIB = self.UIB
        omegaklambda = self.omegaklambda
        cutoff = self.cutoff
        KXs = self.grid.KXs
        KYs = self.grid.KYs
        dJU = self.dJU
        Sigma2 = 0

        for lambda_ in range(1, cutoff):
            for kx in range(Lx):
                for ky in range(Ly):
                    for lambda_1 in range(1, cutoff):
                        for px in range(Lx):
                            for py in range(Ly):
                                W1 = self.vertices.W(kx, ky, lambda_, px, py, lambda_1)
                                W2 = self.vertices.W(px, py, lambda_1, kx, ky, lambda_)
                                eps = self.epsI(KXs[kx] + KXs[px], KYs[ky] + KYs[py])
                                Sigma2 = Sigma2 + pow(UIB, 2) / 2 * abs(
                                    pow(W1 + W2, 2)
                                ) / (
                                    -omegaklambda[lambda_][kx][ky]
                                    - omegaklambda[lambda_1][px][py]
                                    - dJU * eps
                                    + 1j * 0.0001
                                )
        return Sigma2 / (M * M)

    # def selfenergy(self, kx, ky, lambda_, vks, n0):
    #     dkxs = self.grid.dkxs
    #     dkys = self.grid.dkys
    #     cutoff = self.cutoff
    #     L = self.grid.L
    #     UIB = self.UIB
    #     omegaklambda = self.omegaklambda
    #     cutoff = self.cutoff
    #     KXs = self.grid.KXs
    #     KYs = self.grid.KYs
    #     uks = self.vertices.uks
    #     vks = self.vertices.vks
    #     dJU = self.groundstate.dJU
    #     deltan2 = 0
    #     Sigma1, Sigma2 = 0, 0
    #     res = 0
    #     for lambda_ in range(1, cutoff):
    #         for kx in range(L + 1):
    #             for ky in range(L + 1):
    #                 deltan2 = deltan2 + 1.0 / (4 * np.pi * np.pi) * dkxs[kx] * dkys[
    #                     ky
    #                 ] * self.V(kx, ky, lambda_, kx, ky, lambda_)

    #                 # Sigma1 first term of self energy
    #                 Nkres = self.vertices.Nk(kx, ky, lambda_)
    #                 res = np.real(
    #                     pow(UIB / (2 * np.pi), 2)
    #                     * dkxs[kx]
    #                     * dkys[ky]
    #                     * abs(pow(Nkres, 2))
    #                     / (
    #                         -omegaklambda[lambda_, kx, ky]
    #                         - dJU * epsI(KXs[kx], KYs[ky])
    #                         + 1j * 0.0001
    #                     )
    #                 )
    #                 Sigma1 = Sigma1 + res

    #                 # Sigma1 first term of self energy
    #                 for lambda_1 in range(1, cutoff):
    #                     for px in range(L + 1):
    #                         for py in range(L + 1):
    #                             W1 = W(kx, ky, lambda_, px, py, lambda_1, uks, vks, n0)
    #                             W2 = W(px, py, lambda_1, kx, ky, lambda_, uks, vks, n0)
    #                             Sigma2 = Sigma2 + np.real(
    #                                 pow(UIB, 2)
    #                                 / (2 * pow(2 * np.pi, 4))
    #                                 * (dkxs[kx] * dkys[ky] * dkxs[px] * dkys[py])
    #                                 * abs(pow(W1 + W2, 2))
    #                                 / (
    #                                     -omegaklambda[lambda_, kx, ky]
    #                                     - omegaklambda[lambda_1, px, py]
    #                                     - dJU
    #                                     * epsI(KXs[kx] + KXs[px], KYs[ky] + KYs[py])
    #                                     + 1j * 0.0001
    #                                 )
    #                             )

    #         sigma0s = UIB * (n0 + deltan2)

    #         sigma1s = np.real(Sigma1)

    #         sigma2s = np.real(Sigma2)
    #         return sigma0s, sigma1s, sigma2s


# print("sig0 = ", sigma0s[count], " sig1 = ", sigma1s[count], " sig2 = ", sigma2s[count])

# plt.plot(dJUs, sigma0s + sigma1s + sigma2s )
# plt.show()
