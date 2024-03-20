# Class to define Momentum Grid

import numpy as np

class Grid:
    def __init__(self, Lx, Ly):
        self.Lx = Lx
        self.Ly = Ly
        self.dkx = 2 * np.pi / Lx
        self.dky = 2 * np.pi / Ly
        self.KXs = np.linspace(-np.pi + self.dkx / 2, np.pi - self.dkx / 2, self.Lx)
        self.KYs = np.linspace(-np.pi + self.dky / 2, np.pi - self.dky / 2, self.Ly)
        self.M = self.Lx * self.Ly
        self.dkxs = np.ones(Lx) * self.dkx # does this get used anymore?
        self.dkys = np.ones(Ly) * self.dky
        self.dkxs[0] = self.dkxs[Lx - 1] = self.dkx / 2 # check this.
        self.dkys[0] = self.dkys[Lx - 1] = self.dky / 2
