# All parameters of the model are defined here


class params:
    def __init__(self, max_iter, N, dJU, muU, UIB, cutoff):
        self.max_iter = max_iter
        self.N = N
        self.dJU = dJU
        self.muU = muU
        self.UIB = UIB
        self.cutoff = cutoff
