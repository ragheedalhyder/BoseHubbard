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

    # Write a code that calculates muU for a given sigma0

    muU_qcorr = np.zeros(len(dJUs))
    n0 = 1
    dJUmax = (np.sqrt(n0 + 1) - np.sqrt(n0)) ** 2

    # start with large values of dJU and decrease it
    for count in range(len(dJUs)):
        dJU = dJUs[count]
        found = False
        muU = 0.4
        if dJU < dJUmax:
            muU = np.sqrt(2) - 1
        else:
            muU = muU - 0.015
            while (muU < 0.5) and (
                not found
            ):  # change muU until 1 - pert.sigma0(n0) / UIB = 0
                params = Params(N, dJU, muU, UIB, cutoff)
                gs = groundstate(params)
                cns = gs.cns()
                n0 = gs.n0(cns)
                exc = excitations(grid, params, gs, cns)
                uks, vks, omegaklambda = exc.calculate_matrices()
                verts = vertices(grid, gs, uks, vks, cns, n0)
                pert = perturbative(grid, params, verts, omegaklambda)
                eq = pert.sigma0(n0) / UIB
                print(muU, eq)
                if 1 - eq < 0.0001:
                    found = True
                    break
                muU += 0.0001

        muU_qcorr[count] = muU
        print("dJU = ", dJU, "muU = ", muU_qcorr[count])

    plt.plot(dJUs, muU_qcorr, label=r"$\omega_{0}$")
    plt.show()


if __name__ == "__main__":
    main()
