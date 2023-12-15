import numpy as np

class vertices:
    def __init__(self, grid, groundstate, uks, vks, cns, n0):
        self.grid = grid
        self.groundstate = groundstate
        self.uks = uks
        self.vks = vks
        self.cns = cns
        self.n0 = n0

    def U(self, kx, ky, lambda_, qx, qy, lambda1):
        Ures = 0
        uks = self.uks
        n0 = self.n0
        N = self.groundstate.N
        ns = np.arange(0, N, 1)
        Ures = sum((ns - n0 * (1 - (lambda_ == lambda1) * (lambda_ == 0))) * uks[:][kx][ky][lambda_] * uks[:][qx][qy][lambda1])
        return Ures

    def V(self, kx, ky, lambda_, qx, qy, lambda1):
        Vres = 0
        N = self.groundstate.N
        vks = self.vks
        n0 = self.n0
        ns = np.arange(0, N, 1)
        Vres = sum((ns - n0) * vks[:, kx, ky, lambda_] * vks[:, qx, qy, lambda1])
        return Vres

    def W(self, kx, ky, lambda_, qx, qy, lambda1):
        Wres = 0
        N = self.groundstate.N
        uks = self.uks
        vks = self.vks
        n0 = self.n0
        ns = np.arange(0, N, 1)
        Wres = sum((ns - n0) * uks[:, kx, ky, lambda_] * vks[:, qx, qy, lambda1] + (ns - n0) * uks[:, qx, qy, lambda1] * vks[:, kx, ky, lambda_])
        return Wres

    def Nk(self, kx, ky, lambda_):
        Nkres = 0
        N = self.groundstate.N
        uks = self.uks
        vks = self.vks
        cns = self.cns
        ns = np.arange(0, N, 1)
        for i in range(N):
            Nkres = Nkres + ns[i] * cns[i] * (uks[i, kx, ky, lambda_] + vks[i, kx, ky, lambda_])
        return Nkres
