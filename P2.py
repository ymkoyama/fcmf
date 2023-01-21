import mpmath

class phi:
    def __init__(self, r):
        self.r = mpmath.mpf(r)
        self.sr = mpmath.sqrt(self.r)
        self.A = mpmath.ldexp(1 - self.sr, -1)
        self.B = mpmath.ldexp(1 + self.sr, -1)
        self.rho0 = (mpmath.sqrt(1 + self.r) + mpmath.sqrt(2 * self.sr)) / (1 - self.sr)

    @staticmethod
    def name():
        return "P2"

    @staticmethod
    def latexname(r):
        return r"\varphi_{{\mathcal{{P}}_2,{}}}".format(r)

    def value(self, u):
        y = self.A * u + self.B
        return y * y

    def deriv1(self, u):
        y = self.A * u + self.B
        return 2 * self.A * y

    def hatrho(self):
        return self.rho0

