import numpy as np
from scipy.linalg import inv
from scipy import linalg

class t_mat:
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

    def SigmaPolaron(Epol, KXs, KYs, dkxs, dkys, uks, vks, omegaklambda, dJU, n0, UIB, cutoff):
        M = len(KXs)
        N = len(uks)
        dim = M * M * (cutoff + 1)
        eta = 0.0001
        Us = np.zeros((dim, dim))
        Vs = np.zeros((dim, dim))
        Ws = np.zeros((dim, dim))
        PI1100 = PI1121 = PI21 = PI12 = PI2200 = PI22 = inv11 = inv12 = inv22 = T1100 = T1200 = T1121 = T1222 = T2100 = T2200 = T12 = T22 = np.zeros((dim, dim), dtype=complex)
        indk0 = 0
        indq0 = 0
        ZeroIndx1 = True
        ZeroIndx2 = True
        for lambda_ in range(cutoff + 1):
            for lambda1 in range(cutoff + 1):
                for kx1 in range(M):
                    for ky1 in range(M):
                        for qx1 in range(M):
                            for qy1 in range(M):
                                kx = kx1
                                ky = ky1
                                qx = qx1
                                qy = qy1
                                indk = lambda_ * M * M + kx * M + ky
                                indq = lambda1 * M * M + qx * M + qy
                                if lambda_ == 0 and KXs[kx] == 0 and KYs[ky] == 0 and ZeroIndx1 == True:
                                    indk0 = indk
                                    ZeroIndx1 = False
                                if lambda1 == 0 and KXs[qx] == 0 and KYs[qy] == 0 and ZeroIndx1 == True:
                                    indq0 = indq
                                    ZeroIndx2 = False
                                Uelem = UIB * U(kx, ky, lambda_, qx, qy, lambda1, uks, KXs, KYs, n0)
                                Velem = UIB * V(kx, ky, lambda_, qx, qy, lambda1, vks, KXs, KYs, n0)
                                Welem = UIB * (W(kx, ky, lambda_, qx, qy, lambda1, uks, vks, KXs, KYs, n0) + W(qx, qy, lambda1, kx, ky, lambda_, uks, vks, KXs, KYs, n0))
                                epsplus = epsI(KXs[kx] + KXs[qx], KXs[ky] + KXs[qy])
                                epsminus = epsI(KXs[kx] - KXs[qx], KXs[ky] - KXs[qy])
                                if lambda1 == 0:
                                    epsplus = epsI(KXs[kx], KXs[ky])
                                    epsminus = epsI(KXs[kx], KXs[ky])
                                if lambda_ == 0:
                                    epsplus = epsI(KXs[qx], KXs[qy])
                                    epsminus = epsI(-KXs[qx], -KXs[qy])
                                Den11 = Epol - omegaklambda[lambda1][qx][qy] - dJU * epsI(KXs[qx], KXs[qy]) + eta * 1j
                                Den22 = Epol - omegaklambda[lambda1][qx][qy] - dJU * epsI(KXs[qx], KXs[qy]) + eta * 1j
                                Den = Epol - omegaklambda[lambda_][kx][ky] - omegaklambda[lambda1][qx][qy] - dJU * epsplus + eta * 1j
                                Us[indk,indq] = Uelem
                                Vs[indk,indq] = Velem
                                Ws[indk,indq] = Welem
                                coeff = 1 / (4 * np.pi * np.pi) * dkxs[qx] * dkys[qy]
                                PI1100[indk,indq] = coeff * (Uelem / Den11)
                                PI1121[indk,indq] = coeff * (Uelem / Den22)
                                PI21[indk,indq] = coeff * (Welem / Den22)
                                PI2200[indk,indq] = coeff * (Welem / Den22)
                                PI12[indk,indq] = coeff * (Uelem / Den)