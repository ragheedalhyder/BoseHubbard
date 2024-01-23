import numpy as np


class vertices:
    def __init__(self, grid, params, uks, vks, cns, n0):
        self.grid = grid
        self.N = params.N
        self.uks = uks
        self.vks = vks
        self.cns = cns
        self.n0 = n0

    def Nk(self, kx, ky, lambda_):
        Nkres = 0
        N = self.N
        uks = self.uks
        vks = self.vks
        cns = self.cns
        ns = np.arange(0, N, 1)
        Nkres = sum(ns * cns * (uks[:, kx, ky, lambda_] + vks[:, kx, ky, lambda_]))
        # print("Nkres = ", Nkres, "kx = ", kx, "ky = ", ky, "lambda_ = ", lambda_, "uks = ", uks[:,kx,ky,lambda_], "vks = ", vks[:,kx,ky,lambda_])
        return Nkres

    def U(self, kx, ky, lambda_, qx, qy, lambda1):
        Ures = 0
        uks = self.uks
        n0 = self.n0
        N = self.N
        ns = np.arange(0, N, 1)
        Ures = sum(
            (ns - n0 * (1 - (lambda_ == lambda1) * (lambda_ == 0)))
            * uks[:, kx, ky, lambda_]
            * uks[:, qx, qy, lambda1]
        )
        return Ures

    def V(self, kx, ky, lambda_, qx, qy, lambda1):
        Vres = 0
        N = self.N
        vks = self.vks
        n0 = self.n0
        # print("n0 = ", n0)
        ns = np.arange(0, N, 1)
        Vres = sum((ns - n0) * vks[:, kx, ky, lambda_] * vks[:, qx, qy, lambda1])
        return Vres

    def W(self, kx, ky, lambda_, qx, qy, lambda1):
        Wres = 0
        N = self.N
        uks = self.uks
        vks = self.vks
        n0 = self.n0
        ns = np.arange(0, N, 1)
        Wres = sum(
            (ns - n0) * uks[:, kx, ky, lambda_] * vks[:, qx, qy, lambda1]
        )
        return Wres
