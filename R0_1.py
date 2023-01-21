import mpmath

class phi:
    def __init__(self, r):
        self.r = mpmath.mpf(r)
        self.A = 2 * self.r / (1 - self.r)
        self.B = (1 + self.r) / (1 - self.r)
        sr = mpmath.sqrt(self.r)
        self.rhoA = (1 + sr) / (1 - sr)

    @staticmethod
    def name():
        return "R0_1"

    @staticmethod
    def latexname(r):
        return r"\varphi_{{\mathcal{{R}}_{{0,1}},{}}}".format(r)

    def value(self, u):
        return self.A / (self.B - u)

    def deriv1(self, u):
        d = self.B - u
        return self.A / (d * d)

    def hatrho(self):
        return self.rhoA

