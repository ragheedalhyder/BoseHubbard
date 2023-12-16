import numpy as np
import os
import utils
from grid import grid
from params import params
from gs import groundstate
from exc import excitations
from vert import vertices
from pert import perturbative
from plotting import plot2D, plot_cns, plot_omega0
import matplotlib.pyplot as plt


def main():
    args = utils.parse_args()

    output_dir = utils.create_output_dir()

    config = utils.read_config(args.config)

    Lx = config["grid"]["Lx"]
    Ly = config["grid"]["Ly"]
    max_iter = config["physics"]["max_iter_psi"]
    dJUs = np.arange(config["physics"]["dJU"])
    N = config["grid"]["N"]
    muU = eval(config["physics"]["Mu"])
    UIB = config["physics"]["UIB"]

    grid = grid(Lx, Ly)
    params = params(**config["physics"])
    print(params)
    omega0s = np.zeros(len(dJUs))
    omega1s = np.zeros(len(dJUs))
    omega2s = np.zeros(len(dJUs))
    omegas = np.zeros((len(dJUs), 3))
    for count in range(len(dJUs)):
        dJU = dJUs[count]
        gs = groundstate(max_iter, N, dJU, muU)
        cns = gs.cns()

        exc = excitations(grid, gs, cns)
        uks, vks, omegaklambda = exc.calculate_matrices()
        vert = vertices(grid, gs, uks, vks, cns, gs.n0)
        pert = perturbative(grid, gs, vert, UIB)
        omega0s[count] = pert.sigma0()
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
    plt.plot(dJUs, omega0s, label=r"$\omega_{0}$")
    plt.show()


if __name__ == "__main__":
    main()
