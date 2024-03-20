# We calculate the self-energy for a fixed chemical potential and a range of dJU values

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

config_path = "config.yml" 

output_dir = class_utils.create_output_dir()

config = class_utils.read_config(config_path)

Lx = config["grid"]["Lx"]
Ly = config["grid"]["Ly"]

N = config["physics"]["N"]
dJUs = np.linspace(**config["lists"]["dJUs"])
muU = eval(config["physics"]["muU"])

UIB = config["physics"]["UIB"]
cutoff = config["physics"]["cutoff"]


grid = Grid(Lx, Ly)
io = IO()

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


for dJU_ind in range(0, len(dJUs)):
    dJU = dJUs[dJU_ind]
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
    

io.save_to_hdf5_fixed_chemical_potentials(grid, params, dJUs, en_vector, omega0s, omega1s, omega2s, T11, T12, T21, T22, T22_F, SE_SI)