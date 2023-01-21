import mpmath

class phi:
    def __init__(self, r):
        self.r = mpmath.mpf(r)
        self.a1 = mpmath.ldexp(1 - self.r, -1)
        self.a0 = mpmath.ldexp(1 + self.r, -1)
        sr = mpmath.sqrt(self.r)
        self.rho0 = (1 + sr) / (1 - sr)

    @staticmethod
    def name():
        return "P1"

    @staticmethod
    def latexname(r):
        return r"\varphi_{{\mathcal{{P}}_1,{}}}".format(r)

    def value(self, u):
        return self.a1 * u + self.a0

    def deriv1(self, u):
        return self.a1

    def hatrho(self):
        return self.rho0

    def basis(self, n, x):
        s = mpmath.exp(- self.a0 * x) * mpmath.besseli(n, self.a1 * x)
        return s if n % 2 == 0 else -s

