import numpy as np
import os
import utils
from grid import Grid
from gs import groundstate
from exc import excitations
from vert import vertices
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
    UIB = config["physics"]["UIB"]

    grid = Grid(Lx, Ly)

    for count in range(len(dJUs)):
        dJU = dJUs[count]
        gs = groundstate(max_iter, N, dJU, muU)
        cns = gs.cns()
        n0 = gs.n0(cns)

        exc = excitations(grid, gs, cns)
        uks, vks, omegaklambda = exc.calculate_matrices()

        vert = vertices(grid, gs, uks, vks, cns, n0)
        


if __name__ == "__main__":
    main()
