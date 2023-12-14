import numpy as np
import os
import utils
from grid import Grid
from gs import groundstate
from exc import excitations
from plotting import plot2D, plot_cns, plot_omega0
import matplotlib.pyplot as plt


def main():
    args = utils.parse_args()

    output_dir = utils.create_output_dir()

    config = utils.read_config(args.config)

    Lx = config["grid"]["Lx"]
    Ly = config["grid"]["Ly"]
    max_iter = config["max_iter_psi"]
    dJUs = np.arange(**config["physics"]["dJU"])
    N = config["grid"]["N"]
    muU = eval(config["physics"]["Mu"])

    grid = Grid(Lx, Ly)

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

    #     omega0s[count] = omegaklambda[1][5][5]
    #     omega1s[count] = omegaklambda[2][5][5]
    #     omega2s[count] = omegaklambda[3][5][5]
    # omegas[:, 0] = omega0s
    # omegas[:, 1] = omega1s
    # omegas[:, 2] = omega2s
    # plot_omega0(
    #     omegas,
    #     r"$2\delta J/U$",
    #     r"$\omega_{\lambda}(|\vec{k}| = 0)$",
    #     os.path.join(output_dir, "omega0s.png"),
    # )


if __name__ == "__main__":
    main()
