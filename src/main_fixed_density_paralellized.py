# We calculate the self-energy for a fixed density and a range of dJU values

import numpy as np
import class_utils as class_utils
from class_grid import Grid
from class_params import Params
from class_groundstate import groundstate
from class_excitations import excitations
from class_vertices import vertices
from class_perturbation import perturbative
from class_self_energy import Self_Energy
from class_io import IO
from class_plotting import plot2D, plot_cns, plot_omega0
import matplotlib.pyplot as plt
import concurrent.futures

config_path = "config.yml" 

output_dir = class_utils.create_output_dir()

config = class_utils.read_config(config_path)

Lx = config["grid"]["Lx"]
Ly = config["grid"]["Ly"]

N = config["physics"]["N"]
dJUs = np.linspace(**config["lists"]["dJUs"])
dJU = dJUs[0]
muU = eval(config["physics"]["muU"])
UIB = config["physics"]["UIB"]
cutoff = config["physics"]["cutoff"]

grid = Grid(Lx, Ly)
params = Params(N, dJU, muU, UIB, cutoff)

io = IO()
desired_n0 = 1.0
dJUs, muUs = io.read_from_hdf5_fixed_density_qcorr(grid, params, desired_n0)

en_vector = np.linspace(**config["lists"]["energies"])

omega0s = np.zeros(len(dJUs))
omega1s = np.zeros(len(dJUs))
omega2s = np.zeros(len(dJUs))

SE = np.zeros((7, len(en_vector)), dtype=np.complex128)
SEE = np.zeros((len(en_vector), len(dJUs)), dtype=np.complex128)
T11 = np.zeros((len(dJUs), len(en_vector)), dtype=np.complex128)
T12 = np.zeros((len(dJUs), len(en_vector)), dtype=np.complex128)
T21 = np.zeros((len(dJUs), len(en_vector)), dtype=np.complex128)
T22 = np.zeros((len(dJUs), len(en_vector)), dtype=np.complex128)
T22_F = np.zeros((len(dJUs), len(en_vector)), dtype=np.complex128)
SE_SI = np.zeros((len(dJUs), len(en_vector)), dtype=np.complex128)


def process_dJU(dJU_ind):
    dJU = dJUs[dJU_ind]
    muU = muUs[dJU_ind]
    params = Params(N, dJU, muU, UIB, cutoff)
    gs = groundstate(params)
    cns = gs.cns()
    n0 = gs.n0(cns)
    psi0 = gs.psi0(cns)
    exc = excitations(grid, params, gs, cns)
    uks, vks, omegaklambda = exc.calculate_matrices()
    verts = vertices(grid, gs, uks, vks, cns, n0)
    pert = perturbative(grid, params, verts, omegaklambda)
    Pert_Energy = pert.perturbative_energy(n0)
    omega0s[dJU_ind] = Pert_Energy[0]
    omega1s[dJU_ind] = Pert_Energy[1]
    omega2s[dJU_ind] = Pert_Energy[2]

    for epol_ind in range(len(en_vector)):
        Epol = en_vector[epol_ind]
        self_energy = Self_Energy(Epol, grid, params, verts, omegaklambda)
        output = self_energy.calculate_self_energy()
        SE[:, epol_ind] = output

    T11[dJU_ind, :] = SE[1, :]
    T12[dJU_ind, :] = SE[2, :]
    T21[dJU_ind, :] = SE[3, :]
    T22[dJU_ind, :] = SE[4, :]
    T22_F[dJU_ind, :] = SE[5, :]
    SE_SI[dJU_ind, :] = SE[6, :]
    print("dJU = ", dJU, "muU = ", muU, "n0 = ", n0)

# Create a ThreadPoolExecutor
with concurrent.futures.ThreadPoolExecutor() as executor:
    # Use the executor to map the function to the inputs
    executor.map(process_dJU, range(len(dJUs)))
# save results to hdf5 file
io.save_to_hdf5_all_fixed_density(grid, params, desired_n0, dJUs, muUs, en_vector, omega0s, omega1s, omega2s, T11, T12, T21, T22, T22_F, SE_SI)