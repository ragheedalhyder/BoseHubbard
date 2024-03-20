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

    def U_mat(self):
        uks_2d = np.reshape(self.uks, (self.uks.shape[0], -1), order='F') # (7,700)
        # ‘F’ means to read / write the elements using Fortran-like index order, with the first index changing fastest, and the last index changing slowest
        ns = np.arange(self.N)
        U_2d1 = self.n0 * np.matmul(uks_2d.T, uks_2d) 
        U_2d1[:self.grid.M, :self.grid.M ] = 0  
        U_2d2 = np.matmul(ns * uks_2d.T , uks_2d)
        U_2d = U_2d2 - U_2d1
        return U_2d
    
    def V_mat(self):
        vks_2d = np.reshape(self.vks, (self.vks.shape[0], -1), order='F')
        ns = np.arange(self.N)
        V_2d1 = self.n0 * np.matmul(vks_2d.T, vks_2d)
        V_2d2 = np.matmul(ns * vks_2d.T , vks_2d)
        V_2d = V_2d2 - V_2d1
        return V_2d
    
    def W_mat(self):
        uks_2d = np.reshape(self.uks, (self.uks.shape[0], -1), order='F')
        vks_2d = np.reshape(self.vks, (self.vks.shape[0], -1), order='F')
        ns = np.arange(self.N)
        W_2d1 = self.n0 * np.matmul(uks_2d.T, vks_2d)
        W_2d2 = np.matmul(ns * uks_2d.T , vks_2d)
        W_2d = W_2d2 - W_2d1
        return W_2d
