import numpy as np
from scipy.linalg import inv
from scipy import linalg

class Self_Energy:
    def __init__(self, Epol, grid, params, vertices, omegaklambda):
        self.grid = grid
        self.Lx = grid.Lx
        self.Ly = grid.Ly
        self.N = params.N
        self.M = self.Lx * self.Ly
        self.dJU = params.dJU
        self.muU = params.muU
        self.UIB = params.UIB
        self.cutoff = params.cutoff
        self.vertices = vertices
        self.omegaklambda = omegaklambda
        self.Epol = Epol
        self.eta = 0.001
    
    def __str__(self):
        return f"groundstate = {self.groundstate}, UIB = {self.UIB}, cutoff = {self.cutoff}"
    
    def epsI(self, kx, ky):
        return pow(np.sin(kx / 2), 2) + pow(np.sin(ky / 2), 2) #+ 1 # note the 1!  This is to account for the minimum of the tight-binding band.  Important!
    
    def epsI_vec(self):
        kx_grid = np.repeat(self.grid.KXs, self.grid.Lx)
        ky_grid = np.tile(self.grid.KYs, self.grid.Ly)
        epsI_grid = self.epsI(kx_grid, ky_grid)
        epsI_vec = np.tile(epsI_grid, self.N)
        return epsI_vec
    
    def omega_vec(self):
        omegaklambda = self.omegaklambda
        omega_1d = omegaklambda.ravel()
        return omega_1d

    def omega_mat(self):
        omegaklambda = self.omegaklambda
        omega_1d = omegaklambda.ravel()
        omegaklambda_2d = omega_1d[:, np.newaxis] + omega_1d[np.newaxis, :]
        return omegaklambda_2d
    
    def Den1_vec(self):
        dim = self.grid.M * self.N
        omega_vec = self.omega_vec()
        epsI_vec = self.epsI_vec()
        dJU = self.dJU
        eta = self.eta
        Den1 = self.Epol - omega_vec - dJU * epsI_vec + eta * 1j
        Den1 = np.tile(Den1, (dim,1)).T
        return Den1
    
    def eps_grid(self):
        kx_vec =  np.tile(np.repeat(self.grid.KXs, self.grid.Lx), self.N)
        ky_vec = np.tile(self.grid.KYs, self.grid.Ly * self.N)
        
        epsI_grid = self.epsI(kx_vec[:, np.newaxis ] + kx_vec[np.newaxis, : ], ky_vec[:, np.newaxis ] + ky_vec[np.newaxis, : ])
        return epsI_grid
    
    def delete_elements(self, array, indices, axis = None): # checked
        if axis == None:
            array = np.delete(array, indices)
            return array
        else:
            array = np.delete(array, indices, 0)
            array = np.delete(array, indices, 1)
            return array

    def calculate_self_energy(self):
        M = self.grid.M
        UIB = self.UIB
        dJU = self.dJU

        dim = M * self.N
        omega_vec = self.omega_vec()
        epsI_vec = self.epsI_vec()
        dJU = self.dJU
        eta = self.eta
        
        # a bit mysterious this part...
        # 1-w propagator
        Den1 = self.Epol - omega_vec - dJU * epsI_vec + eta * 1j
        Den1 = np.tile(Den1, (dim,1))

        # 2-w propagator
        omega_mat = omega_vec[:, np.newaxis] + omega_vec[np.newaxis, :]
        epsI_grid = self.eps_grid()
        epsI_grid[0, :] = self.epsI_vec()
        epsI_grid[:, 0] = self.epsI_vec().T
        Den2 = self.Epol - omega_mat - dJU * epsI_grid + eta * 1j
        Den2[0, :] = self.Epol - self.omega_vec() - dJU * self.epsI_vec() + eta * 1j # why is there an epsI_vec and a self.epsI_vec
        # Den2[0, :] = self.Epol + eta * 1j
        Den2[:, 0] = Den2[0, :].T # ?


        U_mat = UIB * self.vertices.U_mat()
        V_mat = UIB * self.vertices.V_mat()
        W_mat = UIB * (self.vertices.W_mat() + self.vertices.W_mat().T) # can be optimized further
	
        Den1 = self.delete_elements(Den1, np.arange(1,self.M), axis = 0)
        Den2 = self.delete_elements(Den2, np.arange(1,self.M), axis = 0)
        U_mat = self.delete_elements(U_mat, np.arange(1,self.M), axis = 0)
        V_mat = self.delete_elements(V_mat, np.arange(1,self.M), axis = 0)
        W_mat = self.delete_elements(W_mat, np.arange(1,self.M), axis = 0)

        Pair_Prop_11_12 = np.divide(U_mat, Den1) / M
        Pair_Prop_21_22 = np.divide(W_mat, Den1) / M
        Pair_Prop_12_SE = np.divide(U_mat, Den2) / M
        Pair_Prop_22_SE = np.divide(W_mat, Den2) / M
        

        # eliminate internal condensate lines
	# set the zerot column to zero.  I'm not sure that's right.  
        Pair_Prop_11_12[:, 0] = 0 # These are not actually equal in this block!
        Pair_Prop_21_22[:, 0] = 0
        Pair_Prop_12_SE[:, 0] = 0
        Pair_Prop_22_SE[:, 0] = 0
#        Pair_Prop_11_12[0,:] = 0 # These are not actually equal in this block!
#        Pair_Prop_21_22[0,:] = 0
#        Pair_Prop_12_SE[0,:] = 0
#        Pair_Prop_22_SE[0,:] = 0

        IMat = np.eye(len(Den1))

        inv11 = inv(IMat - Pair_Prop_11_12) # (601 by 601)
        T_11 = inv11 @ U_mat
        T_12 = inv11 @ W_mat
        T_21 = W_mat + Pair_Prop_21_22 @ T_11
        T_22 = V_mat + Pair_Prop_21_22 @ T_12
        
        inv12 = inv(IMat - Pair_Prop_12_SE)
        T_12_SE = inv12 @ W_mat
        T_22_SE = V_mat + 0.5 * Pair_Prop_22_SE @ T_12_SE

        T22_diag = np.diag(T_22_SE)
        Sigma_22 = sum(T22_diag[1:]) / M
 
        indk0 = 0
        indq0 = 0
        sigpol = T_11[indk0, indq0] + T_12[indk0, indq0] + T_21[indk0, indq0] + T_22[indk0, indq0] + Sigma_22
        return np.array([np.real(self.Epol - sigpol), T_11[indk0, indq0], T_12[indk0, indq0], T_21[indk0, indq0], T_22[indk0, indq0], Sigma_22, sigpol], dtype=np.complex128)


# function for calculating hte self energy to O(UIB^2) for debugging purposes
    def calculate_self_energy_perturbative(self):
        M = self.grid.M
        UIB = self.UIB
        dJU = self.dJU

        dim = M * self.N
        omega_vec = self.omega_vec()
        epsI_vec = self.epsI_vec()
        dJU = self.dJU
        eta = self.eta
        
        
        # a bit mysterious this part...
        # 1-w propagator
        Den1 = self.Epol - omega_vec - dJU * epsI_vec + eta * 1j
        Den1 = np.tile(Den1, (dim,1))

        # 2-w propagator
        omega_mat = omega_vec[:, np.newaxis] + omega_vec[np.newaxis, :]
        epsI_grid = self.eps_grid()
        epsI_grid[0, :] = self.epsI_vec()
        epsI_grid[:, 0] = self.epsI_vec().T
        Den2 = self.Epol - omega_mat - dJU * epsI_grid + eta * 1j
        Den2[0, :] = self.Epol - self.omega_vec() - dJU * self.epsI_vec() + eta * 1j # why is there an epsI_vec and a self.epsI_vec
        # Den2[0, :] = self.Epol + eta * 1j
        Den2[:, 0] = Den2[0, :].T # ?
        
        #Den's are (700,700 arrays)
        
 
 

        U_mat = UIB * self.vertices.U_mat()
        V_mat = UIB * self.vertices.V_mat()
        W_mat = UIB * (self.vertices.W_mat() + self.vertices.W_mat().T) # can be optimized further


# some of these deletions are to remove the internal condensate lines 
# some are to remove equal contributions
# ragheed is skeptical of them.  

# In the lambda = lambda' = 0 block all the elements are the same [tbc]
#        print(Den1.shape)
        Den1 = self.delete_elements(Den1, np.arange(1,self.M), axis = 0)
        Den2 = self.delete_elements(Den2, np.arange(1,self.M), axis = 0)
        U_mat = self.delete_elements(U_mat, np.arange(1,self.M), axis = 0)
        V_mat = self.delete_elements(V_mat, np.arange(1,self.M), axis = 0)
        W_mat = self.delete_elements(W_mat, np.arange(1,self.M), axis = 0)
#        print(Den1.shape)
#        exit()

        Pair_Prop_11_12 = np.divide(U_mat, Den1) / M # U over 1p denom
        Pair_Prop_21_22 = np.divide(W_mat, Den1) / M # W over 1p denom
        Pair_Prop_12_SE = np.divide(U_mat, Den2) / M # U over 2p denom 
        Pair_Prop_22_SE = np.divide(W_mat, Den2) / M # W over 2p denom
        
        # eliminate internal condensate lines
	# set the zerot column to zero.  I'm not sure that's right.  
        Pair_Prop_11_12[:, 0] = 0 # These are not actually equal in this block!
        Pair_Prop_21_22[:, 0] = 0
        Pair_Prop_12_SE[:, 0] = 0
        Pair_Prop_22_SE[:, 0] = 0
        
#        Pair_Prop_11_12[:, 0:M] = 0 # These are not actually equal in this block!
#        Pair_Prop_21_22[:, 0:M] = 0
#        Pair_Prop_12_SE[:, 0:M] = 0
#        Pair_Prop_22_SE[:, 0:M] = 0


#        Pair_Prop_11_12[0,:] = 0 # These are not actually equal in this block!
#        Pair_Prop_21_22[0,:] = 0
#        Pair_Prop_12_SE[0,:] = 0
#        Pair_Prop_22_SE[0,:] = 0

        IMat = np.eye(len(Den1))

        #inv11 = inv(IMat - Pair_Prop_11_12) # (601 by 601)
        T_11 = U_mat + Pair_Prop_11_12 @ U_mat
        T_12 = W_mat + Pair_Prop_11_12 @ W_mat
        T_21 = W_mat + Pair_Prop_21_22 @ U_mat #T_11
        T_22 = V_mat + Pair_Prop_21_22 @ W_mat#T_12
        
        #inv12 = inv(IMat - Pair_Prop_12_SE)
        T_12_SE = W_mat +  Pair_Prop_12_SE @ W_mat
        T_22_SE = V_mat + 0.5 * Pair_Prop_22_SE @ W_mat

        T22_diag = np.diag(T_22_SE)
        Sigma_22 = sum(T22_diag[1:]) / M
#        Sigma_22 = sum(T22_diag[M:]) / M
 
        indk0 = 0
        indq0 = 0
        sigpol = T_11[indk0, indq0] + T_12[indk0, indq0] + T_21[indk0, indq0] + T_22[indk0, indq0] + Sigma_22
        return np.array([np.real(self.Epol - sigpol), T_11[indk0, indq0], T_12[indk0, indq0], T_21[indk0, indq0], T_22[indk0, indq0], Sigma_22, sigpol], dtype=np.complex128)




































