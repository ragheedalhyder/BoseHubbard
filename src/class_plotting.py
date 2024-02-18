from matplotlib import pyplot as plt


def plot2D(data, xlabel, ylabel, output_file, title=None):
    plt.plot(data)
    plt.grid()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if title:
        plt.title(title)
    plt.savefig(output_file)


def plot_cns(cnslist, xlabel, ylabel, output_file, title=None):
    for cns in cnslist:
        plt.plot(cns)
    plt.grid()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if title:
        plt.title(title)
    plt.savefig(output_file)


def plot_omega0(omegaslist, xlabel, ylabel, output_file=None, title=None, show=None):
    for omegas in omegaslist:
        plt.plot(omegas)
    plt.grid()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if title:
        plt.title(title)
    if output_file:
        plt.savefig(output_file)
    if show:
        plt.show()
