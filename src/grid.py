# Class to define Momentum Grid

import numpy as np


class Grid:
    def __init__(self, Lx, Ly):
        self.Lx = Lx
        self.Ly = Ly
        self.dkx = 2 * np.pi / Lx
        self.dky = 2 * np.pi / Ly
        self.KXs = np.linspace(-np.pi + self.dkx, np.pi - self.dkx, self.Lx)
        self.KYs = np.linspace(-np.pi + self.dky, np.pi - self.dky, self.Ly)
        self.M = self.Lx * self.Ly
        self.dkxs = np.ones(Lx) * self.dkx
        self.dkys = np.ones(Ly) * self.dky
        self.dkxs[0] = self.dkxs[Lx - 1] = self.dkx / 2
        self.dkys[0] = self.dkys[Lx - 1] = self.dky / 2
