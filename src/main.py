import numpy as np
import os
import utils
from grid import Grid
from cns import psi
from plotting import plot2D, plot_cns


def main():
    args = utils.parse_args()

    output_dir = utils.create_output_dir()

    config = utils.read_config(args.config)

    grid = Grid(config["grid"]["Lx"], config["grid"]["Ly"])

    dJUs = np.arange(**config["physics"]["dJU"])
    print(dJUs)
    cnslist = []
    for dJU in dJUs:
        cns = psi(config["max_iter_psi"], config["grid"]["N"], dJU, eval(config["physics"]["Mu"]))
        cnslist.append(cns)

    plot_cns(cnslist, r'$n$', r'$\bar{c}_n$', os.path.join(output_dir, "csn.png") )


if __name__ == "__main__":
    main()
