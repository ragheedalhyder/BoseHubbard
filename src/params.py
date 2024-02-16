# All parameters of the model are defined here


class Params:
    def __init__(self, N, dJU, muU, UIB, cutoff, Epol = None, **kwargs):
        self.N = N
        self.dJU = dJU
        self.muU = muU
        self.UIB = UIB
        self.cutoff = cutoff
        if (Epol is None):
            self.Epol = []
        else:
            self.Epol = Epol
