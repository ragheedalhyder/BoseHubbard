# Perturbative diagrams
import numpy as np

class perturbative:
    def __init__(self, grid, groundstate, vertices, omegaklambda, UIB, cutoff):
        self.grid = grid
        self.groundstate = groundstate
        self.vertices = vertices
        self.UIB = UIB
        self.cutoff = cutoff
        self.omegaklambda = omegaklambda

    def __str__(self):
        return f"groundstate = {self.groundstate}, UIB = {self.UIB}, cutoff = {self.cutoff}"
    
    def __repr__(self):
        return f"groundstate = {self.groundstate}, UIB = {self.UIB}, cutoff = {self.cutoff}"
    
    # Sigma0 zeroth term of self energy
    def sigma0(self, kx, ky, lambda_, n0):
        deltan2 = 0
        N = self.groundstate.N
        dkxs = self.grid.dkxs
        dkys = self.grid.dkys
        cutoff = self.cutoff
        L = self.grid.L
        UIB = self.UIB
        cutoff = self.cutoff

        for lambda_ in range(1, cutoff):
            for kx in range(L + 1):
                for ky in range(L + 1):
                    V_vert = self.vertices.V(kx, ky, lambda_, kx, ky, lambda_)
                    deltan2 = deltan2 + 1.0 / (4 * np.pi * np.pi) * dkxs[kx] * dkys[ky] * V_vert
        return UIB * ( n0 + deltan2 )
    

    def sigma1(self, kx, ky, lambda_, n0):
        deltan2 = 0
        cns = self.vertices.cns
        dkxs = self.grid.dkxs
        dkys = self.grid.dkys
        cutoff = self.cutoff
        L = self.grid.L
        UIB = self.UIB
        omegaklambda = self.omegaklambda
        cutoff = self.cutoff
        dkxs = self.grid.dkxs
        dkys = self.grid.dkys
        KXs = self.grid.KXs
        KYs = self.grid.KYs
        dJU = self.groundstate.dJU

        for lambda_ in range(1, cutoff):
            for kx in range(L + 1):
                for ky in range(L + 1):
                    Nkres = self.vertices.Nk(kx, ky, lambda_)
                    eps = self.groundstate.epsI(KXs[kx], KYs[ky])
                    res = np.real( pow( UIB / ( 2 * np.pi ) , 2 ) * dkxs[kx] * dkys[ky] * abs(pow(Nkres, 2)) / ( - omegaklambda[lambda_ , kx , ky] - dJU * eps + 1j * 0.0001))
                    Sigma1 = Sigma1 + res
        return UIB * (n0 + deltan2)
    
    def sigma2(self, kx, ky, lambda_, vks, n0):
        dkxs = self.grid.dkxs
        dkys = self.grid.dkys
        cutoff = self.cutoff
        L = self.grid.L
        UIB = self.UIB
        omegaklambda = self.omegaklambda
        cutoff = self.cutoff
        KXs = self.grid.KXs
        KYs = self.grid.KYs
        uks = self.vertices.uks
        vks = self.vertices.vks
        dJU = self.groundstate.dJU

        for lambda_ in range(1, cutoff):
            for kx in range(L + 1):
                for ky in range(L + 1):
                    for lambda_1 in range(1, cutoff):
                        for px in range(L + 1):
                            for py in range(L + 1):
                                W1 = self.W(kx, ky, lambda_, px, py, lambda_1)
                                W2 = self.W(px, py, lambda_1, kx, ky, lambda_)
                                eps = self.groundstate.epsI(KXs[kx] + KXs[px], KYs[ky] + KYs[py])
                                Sigma2 = Sigma2 + np.real(pow(UIB, 2) / (2 * pow(2 * np.pi,4) ) * (dkxs[kx] * dkys[ky] * dkxs[px] * dkys[py]) * abs(pow(W1 + W2, 2)) / (-omegaklambda[lambda_ , kx , ky] - omegaklambda[lambda_1 , px , py] - dJU * eps + 1j * 0.0001))
    
    def selfenergy(self, kx, ky, lambda_, vks, n0):
        dkxs = self.grid.dkxs
        dkys = self.grid.dkys
        cutoff = self.cutoff
        L = self.grid.L
        UIB = self.UIB
        omegaklambda = self.omegaklambda
        cutoff = self.cutoff
        KXs = self.grid.KXs
        KYs = self.grid.KYs
        uks = self.vertices.uks
        vks = self.vertices.vks
        dJU = self.groundstate.dJU
        deltan2 = 0
        Sigma1, Sigma2 = 0 , 0
        res = 0
        for lambda_ in range(1, cutoff):
            for kx in range(L + 1):
                for ky in range(L + 1):
                    deltan2 = deltan2 + 1.0 / (4 * np.pi * np.pi) * dkxs[kx] * dkys[ky] * V(kx, ky, lambda_, kx, ky, lambda_, vks, n0)

                    # Sigma1 first term of self energy
                    Nkres = self.vertices.Nk(kx, ky, lambda_)
                    res = np.real( pow( UIB / ( 2 * np.pi ) , 2 ) * dkxs[kx] * dkys[ky] * abs(pow(Nkres, 2)) / ( - omegaklambda[lambda_ , kx , ky] - dJU * epsI(KXs[kx], KYs[ky]) + 1j * 0.0001))
                    Sigma1 = Sigma1 + res

                    # Sigma1 first term of self energy
                    for lambda_1 in range(1, cutoff):
                        for px in range(L + 1):
                            for py in range(L + 1):
                                W1 = W(kx, ky, lambda_, px, py, lambda_1, uks, vks, n0)
                                W2 = W(px, py, lambda_1, kx, ky, lambda_, uks, vks, n0)
                                Sigma2 = Sigma2 + np.real(pow(UIB, 2) / (2 * pow(2 * np.pi,4) ) * (dkxs[kx] * dkys[ky] * dkxs[px] * dkys[py]) * abs(pow(W1 + W2, 2)) / (-omegaklambda[lambda_ , kx , ky] - omegaklambda[lambda_1 , px , py] - dJU * epsI(KXs[kx] + KXs[px], KYs[ky] + KYs[py]) + 1j * 0.0001))
    

sigma0s[count] = UIB * (n0 + deltan2)

sigma1s[count] = np.real(Sigma1)

sigma2s[count] = np.real(Sigma2)

# print("sig0 = ", sigma0s[count], " sig1 = ", sigma1s[count], " sig2 = ", sigma2s[count])

# plt.plot(dJUs, sigma0s + sigma1s + sigma2s )
# plt.show()