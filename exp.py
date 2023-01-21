import mpmath

class phi:
    def __init__(self, r):
        self.r = mpmath.mpf(r)
        self.sr = mpmath.sqrt(self.r)
        self.C = mpmath.ldexp(mpmath.log(1 / self.r), -1)
        s = mpmath.pi / mpmath.log(1 / self.r)
        self.rho0 = s + mpmath.sqrt(s * s + 1)

    @staticmethod
    def name():
        return "exp"

    @staticmethod
    def latexname(r):
        return r"\varphi_{{\exp,{}}}".format(r)

    def value(self, u):
        return self.sr * mpmath.exp(self.C * u)

    def deriv1(self, u):
        return self.C * self.value(u)

    def hatrho(self):
        return self.rho0

