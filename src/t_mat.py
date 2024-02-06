import numpy as np
from scipy.linalg import inv
from scipy import linalg

class t_mat:
    def __init__(self, Epol, grid, params, vertices, omegaklambda):
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
        self.Epol = Epol
        self.eta = 0.0001
    
    def __str__(self):
        return f"groundstate = {self.groundstate}, UIB = {self.UIB}, cutoff = {self.cutoff}"
    
    def epsI(self, kx, ky):
        return pow(np.sin(kx / 2), 2) + pow(np.sin(ky / 2), 2)
    
    def epsI_vec(self):
        kx_grid, ky_grid = np.meshgrid(self.grid.KXs, self.grid.KYs)
        epsI_grid = self.epsI(kx_grid, ky_grid)
        epsI_vec = np.tile(epsI_grid.ravel(), self.N)
        return epsI_vec
    
    def omega_vec(self):
        omegaklambda = self.omegaklambda
        omega_1d = omegaklambda.ravel()
        return omega_1d

    def omega_mat(self):
        omegaklambda = self.omegaklambda
        omega_1d = omegaklambda.ravel()
        omegaklambda_2d = omega_1d[:, np.newaxis] - omega_1d[np.newaxis, :]
        return omegaklambda_2d
    
    def Den1_vec(self):
        omega_vec = self.omega_vec()
        epsI_vec = self.epsI_vec()
        dJU = self.dJU
        eta = self.eta
        Den1 = self.Epol - omega_vec - dJU * epsI_vec + eta * 1j
        return Den1
    
    def eps_grid(self):
        kx_grid, qx_grid = np.meshgrid(np.tile(np.repeat(self.grid.KXs, self.grid.Lx), self.N), np.tile(np.repeat(self.grid.KXs, self.grid.Lx), self.N))
        ky_grid, qy_grid = np.meshgrid(np.tile(self.grid.KXs, self.grid.Lx * self.N), np.tile(self.grid.KXs, self.Lx * N))
        # Apply the epsI function to each pair of kx and ky values
        epsI_grid = self.epsI(kx_grid + qx_grid, ky_grid + qy_grid)
        return epsI_grid

    def Self_Energy(self):
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

        dim = M * M * (cutoff + 1)
        eta = 0.0001
        omega_vec = self.omega_vec()
        epsI_vec = self.epsI_vec()
        dJU = self.dJU
        eta = self.eta
        Den1 = self.Epol - omega_vec - dJU * epsI_vec + eta * 1j
        Den1 = Den1[:, np.newaxis]
        omega_mat = omega_vec[:, np.newaxis] + omega_vec[np.newaxis, :]
        epsI_grid = self.eps_grid()
        Den2 = self.Epol - omega_mat - dJU * epsI_grid + eta * 1j

        U_mat = UIB * self.vertices.U_mat()
        V_mat = UIB * self.vertices.V_mat()
        W_mat = UIB * self.vertices.W_mat()

        PI1100 = np.divide(U_mat, Den1)
        PI1121[indk, indq] = coeff * (Uelem / Den22)
        PI21[indk, indq] = coeff * (Welem / Den22)
        PI2200[indk, indq] = coeff * (Welem / Den22)
        PI12[indk, indq] = coeff * (Uelem / Den)
        PI22[indk, indq] = coeff * (Welem / Den)

        PI1100 = PI1121 = PI21 = PI12 = PI2200 = PI22 = inv11 = inv12 = inv22 = T1100 = T1200 = T1121 = T1222 = T2100 = T2200 = T12 = T22 = np.zeros((dim, dim), dtype=complex)
        indk0 = 0
        indq0 = 0
        ZeroIndx1 = True
        ZeroIndx2 = True

        for lambda_ in range(cutoff):
            for lambda1 in range(cutoff):
                for kx1 in range(Lx):
                    for ky1 in range(Ly):
                        for qx1 in range(Lx):
                            for qy1 in range(Ly):
                                kx = kx1
                                ky = ky1
                                qx = qx1
                                qy = qy1
                                indk = lambda_ * Lx * Lx + kx * Ly + ky
                                indq = lambda1 * Ly * Ly + qx * Ly + qy
                                if lambda_ == 0 and KXs[kx] == 0 and KYs[ky] == 0 and ZeroIndx1 == True:
                                    indk0 = indk
                                    ZeroIndx1 = False
                                if lambda1 == 0 and KXs[qx] == 0 and KYs[qy] == 0 and ZeroIndx1 == True:
                                    indq0 = indq
                                    ZeroIndx2 = False
                                Uelem = UIB * self.vertices.U(kx, ky, lambda_, qx, qy, lambda1)
                                Velem = UIB * self.vertices.V(kx, ky, lambda_, qx, qy, lambda1)
                                Welem = UIB * (self.vertices.W(kx, ky, lambda_, qx, qy, lambda1) + self.vertices.W(qx, qy, lambda1, kx, ky, lambda_))
                                epsplus = self.epsI(KXs[kx] + KXs[qx], KXs[ky] + KXs[qy])
                                epsminus = self.epsI(KXs[kx] - KXs[qx], KXs[ky] - KXs[qy])

                                if lambda1 == 0:
                                    epsplus = self.epsI(KXs[kx], KXs[ky])
                                    epsminus = self.epsI(KXs[kx], KXs[ky])
                                if lambda_ == 0:
                                    epsplus = self.epsI(KXs[qx], KXs[qy])
                                    epsminus = self.epsI(-KXs[qx], -KXs[qy])
                                Den11 = self.Epol - omegaklambda[lambda1, qx, qy] - dJU * self.epsI(KXs[qx], KYs[qy]) + eta * 1j
                                Den22 = self.Epol - omegaklambda[lambda1, qx, qy] - dJU * self.epsI(KXs[qx], KYs[qy]) + eta * 1j
                                Den = self.Epol - omegaklambda[lambda_, kx, ky] - omegaklambda[lambda1, qx, qy] - dJU * epsplus + eta * 1j
                                Us[indk, indq] = Uelem
                                Vs[indk, indq] = Velem
                                Ws[indk, indq] = Welem
                                coeff = 1 / (4 * np.pi * np.pi) * dkxs[qx] * dkys[qy]
                                PI1100[indk, indq] = coeff * (Uelem / Den11)
                                PI1121[indk, indq] = coeff * (Uelem / Den22)
                                PI21[indk, indq] = coeff * (Welem / Den22)
                                PI2200[indk, indq] = coeff * (Welem / Den22)
                                PI12[indk, indq] = coeff * (Uelem / Den)
                                PI22[indk, indq] = coeff * (Welem / Den)
                                if lambda1 == 0:
                                    PI1100[indk, indq] = 0
                                    PI1121[indk, indq] = 0
                                    PI21[indk, indq] = 0
                                    PI12[indk, indq] = 0
                                    PI2200[indk, indq] = 0
                                    PI22[indk, indq] = 0
        IMat = np.eye(dim)
        inv11 = inv(IMat - PI1100)
        inv22 = inv(IMat - PI1121)
        T1100 = inv11 @ Us
        T1200 = inv11 @ Ws
        T1121 = inv22 @ Us
        T2100 = Ws + PI21 @ T1121
        T1222 = inv22 @ Ws
        T2200 = Vs + PI21 @ T1222
        inv12 = inv(IMat - PI12)
        T12 = inv12 @ Ws
        T22 = Vs + 0.5 * PI22 @ T12
        T22diag = np.diag(T22)
        for lambda_ in range(1, cutoff + 1):
            for kx in range(M):
                for ky in range(M):
                    indk = lambda_ * M * M + kx * M + ky
                    Sigma22 += 1.0 / (4 * np.pi * np.pi) * dkxs[kx] * dkxs[ky] * T22diag[indk]
        sigpol = T1100[indk0, indq0] + T1200[indk0, indq0] + T2100[indk0, indq0] + T2200[indk0, indq0] + Sigma22
        return np.array([np.real(Epol - sigpol), T1100[indk0, indq0], T12[indk0, indq0], T2100[indk0, indq0], T2200[indk0, indq0], Sigma22, sigpol], dtype=np.complex128)

    def SigmaPolaron(self):
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
        for lambda_ in range(cutoff):
            for lambda1 in range(cutoff):
                for kx1 in range(Lx):
                    for ky1 in range(Ly):
                        for qx1 in range(Lx):
                            for qy1 in range(Ly):
                                kx = kx1
                                ky = ky1
                                qx = qx1
                                qy = qy1
                                indk = lambda_ * Lx * Lx + kx * Ly + ky
                                indq = lambda1 * Ly * Ly + qx * Ly + qy
                                if lambda_ == 0 and KXs[kx] == 0 and KYs[ky] == 0 and ZeroIndx1 == True:
                                    indk0 = indk
                                    ZeroIndx1 = False
                                if lambda1 == 0 and KXs[qx] == 0 and KYs[qy] == 0 and ZeroIndx1 == True:
                                    indq0 = indq
                                    ZeroIndx2 = False
                                Uelem = UIB * self.vertices.U(kx, ky, lambda_, qx, qy, lambda1)
                                Velem = UIB * self.vertices.V(kx, ky, lambda_, qx, qy, lambda1)
                                Welem = UIB * (self.vertices.W(kx, ky, lambda_, qx, qy, lambda1) + self.vertices.W(qx, qy, lambda1, kx, ky, lambda_))
                                epsplus = self.epsI(KXs[kx] + KXs[qx], KXs[ky] + KXs[qy])
                                epsminus = self.epsI(KXs[kx] - KXs[qx], KXs[ky] - KXs[qy])

                                if lambda1 == 0:
                                    epsplus = self.epsI(KXs[kx], KXs[ky])
                                    epsminus = self.epsI(KXs[kx], KXs[ky])
                                if lambda_ == 0:
                                    epsplus = self.epsI(KXs[qx], KXs[qy])
                                    epsminus = self.epsI(-KXs[qx], -KXs[qy])
                                Den11 = self.Epol - omegaklambda[lambda1, qx, qy] - dJU * self.epsI(KXs[qx], KYs[qy]) + eta * 1j
                                Den22 = self.Epol - omegaklambda[lambda1, qx, qy] - dJU * self.epsI(KXs[qx], KYs[qy]) + eta * 1j
                                Den = self.Epol - omegaklambda[lambda_, kx, ky] - omegaklambda[lambda1, qx, qy] - dJU * epsplus + eta * 1j
                                Us[indk, indq] = Uelem
                                Vs[indk, indq] = Velem
                                Ws[indk, indq] = Welem
                                coeff = 1 / (4 * np.pi * np.pi) * dkxs[qx] * dkys[qy]
                                PI1100[indk, indq] = coeff * (Uelem / Den11)
                                PI1121[indk, indq] = coeff * (Uelem / Den22)
                                PI21[indk, indq] = coeff * (Welem / Den22)
                                PI2200[indk, indq] = coeff * (Welem / Den22)
                                PI12[indk, indq] = coeff * (Uelem / Den)
                                PI22[indk, indq] = coeff * (Welem / Den)
                                if lambda1 == 0:
                                    PI1100[indk, indq] = 0
                                    PI1121[indk, indq] = 0
                                    PI21[indk, indq] = 0
                                    PI12[indk, indq] = 0
                                    PI2200[indk, indq] = 0
                                    PI22[indk, indq] = 0
        IMat = np.eye(dim)
        inv11 = inv(IMat - PI1100)
        inv22 = inv(IMat - PI1121)
        T1100 = inv11 @ Us
        T1200 = inv11 @ Ws
        T1121 = inv22 @ Us
        T2100 = Ws + PI21 @ T1121
        T1222 = inv22 @ Ws
        T2200 = Vs + PI21 @ T1222
        inv12 = inv(IMat - PI12)
        T12 = inv12 @ Ws
        T22 = Vs + 0.5 * PI22 @ T12
        T22diag = np.diag(T22)
        for lambda_ in range(1, cutoff + 1):
            for kx in range(M):
                for ky in range(M):
                    indk = lambda_ * M * M + kx * M + ky
                    Sigma22 += 1.0 / (4 * np.pi * np.pi) * dkxs[kx] * dkxs[ky] * T22diag[indk]
        sigpol = T1100[indk0, indq0] + T1200[indk0, indq0] + T2100[indk0, indq0] + T2200[indk0, indq0] + Sigma22
        return np.array([np.real(Epol - sigpol), T1100[indk0, indq0], T12[indk0, indq0], T2100[indk0, indq0], T2200[indk0, indq0], Sigma22, sigpol], dtype=np.complex128)