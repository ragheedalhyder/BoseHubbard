# All parameters of the model are defined here


class Params:
    def __init__(self, N, dJU, muU, UIB, cutoff, **kwargs):
        self.N = N
        self.dJU = dJU
        self.muU = muU
        self.UIB = UIB
        self.cutoff = cutoff
