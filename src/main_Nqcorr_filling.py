# We calculate fixed density lines for a range of dJU values

import numpy as np
import class_utils as utils
from class_grid import Grid
from class_params import Params
from class_groundstate import groundstate
from class_excitations import excitations
from class_vertices import vertices
from class_perturbation import perturbative
# from class_plotting import plot2D, plot_cns, plot_omega0
# from matplotlib import pyplot as plt
from class_io import IO

config_path = "config.yml" 

output_dir = utils.create_output_dir()

config = utils.read_config(config_path)

Lx = config["grid"]["Lx"]
Ly = config["grid"]["Ly"]

N = 12
dJUs = np.linspace(**config["lists"]["dJUs"])
dJU = dJUs[0]
muU = eval(config["physics"]["muU"])
UIB = config["physics"]["UIB"]
cutoff = config["physics"]["cutoff"]

grid = Grid(Lx, Ly)
params = Params(N, dJU, muU, UIB, cutoff)
io = IO()

desired_n0 = 1. # desired density line
muU_Nqcorr = np.zeros(len(dJUs))
dJUmax = (np.sqrt(desired_n0 + 1) - np.sqrt(desired_n0))**2

def n0_corr(muU, dJU, UIB, cutoff, desired_n0):
    params = Params(N, dJU, muU, UIB, cutoff)
    gs = groundstate(params)
    cns = gs.cns()
    n0 = gs.n0(cns)
    eq = n0
    return eq - desired_n0

def findroot(mu_start, mu_step, dJU, UIB, cutoff, desired_n0):
    n_1 = n0_corr(mu_start, dJU, UIB, cutoff, desired_n0)
    n_2 = n0_corr(mu_start + mu_step, dJU, UIB, cutoff, desired_n0)
    
    if(n_1 * n_2 > 0):
        return findroot(mu_start + mu_step, mu_step, dJU, UIB, cutoff, desired_n0)
    else:
        n_mid = n0_corr(mu_start + mu_step / 2, dJU, UIB, cutoff, desired_n0)
        if(np.abs(n_mid) < 1e-6):
            return mu_start + mu_step / 2
        elif(n_mid * n_1 < 0):
            return findroot(mu_start, mu_step / 2, dJU, UIB, cutoff, desired_n0)
        else:
            return findroot(mu_start + mu_step / 2, mu_step / 2, dJU, UIB, cutoff, desired_n0)


# start with large values of dJU and decrease it 
for count in range(len(dJUs)):
    if(dJUs[count] < dJUmax):
        muU_Nqcorr[count] = np.sqrt(2) - 1
    else:
        dJU = dJUs[count]
        found = False
        mu_start = muU_Nqcorr[count - 1] - 0.3
        mu_step = 0.05
        muU_Nqcorr[count] = findroot(mu_start, mu_step, dJU, UIB, cutoff, desired_n0)
        # print("dJU = ", dJU, "muU = ", muU_qcorr[count])

io.save_to_hdf5_fixed_density_Nqcorr(grid, params, desired_n0, dJUs, muU_Nqcorr)

