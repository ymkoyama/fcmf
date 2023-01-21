import mpmath

class CMF:
    def __init__(self, a, b, param = "1"):
        self.a = mpmath.mpf(a)
        self.b = mpmath.mpf(b)
        self.eta = mpmath.mpf(param)
        self.logba = mpmath.log(self.b / self.a)
        self.etam1 = self.eta - 1
        self.rgamma = mpmath.rgamma(self.eta)

    def w(self, t):
        return self.rgamma * mpmath.power(t, self.etam1)

    def derivn(self, x, n):
        ne = n + self.eta
        fn = self.rgamma if n % 2 == 0 else - self.rgamma
        if x == 0:
            fn *= mpmath.power(self.a, ne) * mpmath.expm1(ne * self.logba)
            fn /= ne
        else:
            fn *= mpmath.gammainc(ne, self.a * x, self.b * x)
            fn /= mpmath.power(x, ne)
        return fn

    def value(self, x):
        return self.derivn(x, 0)

    def deriv1(self, x):
        return self.derivn(x, 1)

    def deriv2(self, x):
        return self.derivn(x, 2)

