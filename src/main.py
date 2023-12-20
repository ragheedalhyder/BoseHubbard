import numpy as np
import os
import utils
from grid import Grid
from params import Params
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

    N = config["physics"]["N"]
    dJUs = np.arange(**config["lists"]["dJUs"])
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
        # print(params.muU)
        gs = groundstate(params)
        cns = gs.cns()
        n0 = gs.n0(cns)
        # print(n0)
        exc = excitations(grid, params, gs, cns)
        uks, vks, omegaklambda = exc.calculate_matrices()
        verts = vertices(grid, gs, uks, vks, cns, n0)
        pert = perturbative(grid, params, verts, omegaklambda)
        omega0s[count] = pert.sigma0(n0)
        omega1s[count] = pert.sigma1()
        omega2s[count] = pert.sigma2()
        print(omega0s[count], omega1s[count], omega2s[count])
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
