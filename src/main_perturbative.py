import numpy as np
import os
import Codes.Python.src.class_utils as class_utils
from Codes.Python.src.class_grid import Grid
from Codes.Python.src.class_params import Params
from Codes.Python.src.class_groundstate import groundstate
from Codes.Python.src.class_excitations import excitations
from Codes.Python.src.class_vertices import vertices
from Codes.Python.src.class_perturbation import perturbative
from Codes.Python.src.class_plotting import plot2D, plot_cns, plot_omega0
import matplotlib.pyplot as plt

def main():
    args = class_utils.parse_args()

    output_dir = class_utils.create_output_dir()

    config = class_utils.read_config(args.config)

    Lx = config["grid"]["Lx"]
    Ly = config["grid"]["Ly"]

    N = config["physics"]["N"]
    dJUs = np.linspace(**config["lists"]["dJUs"])
    muU = eval(config["physics"]["muU"])
    UIB = config["physics"]["UIB"]
    cutoff = config["physics"]["cutoff"]

    grid = Grid(Lx, Ly)

    omega0s = np.zeros(len(dJUs))
    omega1s = np.zeros(len(dJUs))
    omega2s = np.zeros(len(dJUs))

    omegas = np.zeros((len(dJUs), 3))
    for count in range(len(dJUs)):
        dJU = dJUs[count]
        params = Params(N, dJU, muU, UIB, cutoff)
        gs = groundstate(params)
        cns = gs.cns()
        n0 = gs.n0(cns)
        exc = excitations(grid, params, gs, cns)
        uks, vks, omegaklambda = exc.calculate_matrices()
        verts = vertices(grid, gs, uks, vks, cns, n0)
        pert = perturbative(grid, params, verts, omegaklambda)
        Pert_Energy = pert.perturbative_energy(n0)
        omega0s[count] = Pert_Energy[0]
        omega1s[count] = Pert_Energy[1]
        omega2s[count] = Pert_Energy[2]
        # print(omega0s[count], omega1s[count], omega2s[count])
    #     omega0s[count] = omegaklambda[1][5][5]
    #     omega1s[count] = omegaklambda[2][5][5]
    #     omega2s[count] = omegaklambda[3][5][5]
    # omegas[:, 0] = omega0s
    # omegas[:, 1] = omega1s
    # omegas[:, 2] = omega2s
    # plot_omega0(
    #     omegas,
    #     r"$2\delta J/U$",
    #     r"$\omega_{\lambda}(|\vec{k}| = 0)$", show=True
    # )
    plt.plot(dJUs, omega0s + omega1s + omega2s, label=r"$\omega_{0}$")
    plt.show()
    # plt.plot(dJUs, omega1s, label=r"$\omega_{0}$")
    # plt.show()


if __name__ == "__main__":
    main()
