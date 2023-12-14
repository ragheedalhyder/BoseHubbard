import numpy as np
import os
import utils
from grid import Grid
from gs import groundstate
from plotting import plot2D, plot_cns
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

    gs = groundstate(max_iter, N, dJUs[10], muU)

    cns = gs.cns()
    omega0U = gs.omega0U(cns)
    psi0 = gs.psi0(cns)
    n0 = gs.n0(cns)

    plt.plot(cns)
    plt.xlim(0, 70)
    plt.show()
    # psi0 = psi0(cns, N)
    # n0 = n0(cns, N)
    # omega0U = omega0(cns, N, dJU, muU, psi0)

    # n0s[count] = n0
    # omega0U = - 2 * dJU * psi0 * psi0 + omega0U
    # for dJU in dJUs:
    #     cns = psi(config["max_iter_psi"], config["grid"]["N"], 0.5, eval(config["physics"]["Mu"]))
    #     cnslist.append(cns)

    # plot_cns(cnslist, r'$n$', r'$\bar{c}_n$', os.path.join(output_dir, "csn.png") )


if __name__ == "__main__":
    main()
