import mpmath

#return q, pi/(2K(k)), k
def q_pi_2K_k(r):
    a = mpmath.mpf(1)
    g = r
    k2 = 1 - r * r
    k = mpmath.sqrt(k2)
    q = k2 / mpmath.ldexp(r, 4)
    s = mpmath.mpf(3)

    for i in range(10000):
        ai = a
        gi = g
        a = mpmath.ldexp(ai + gi, -1)
        g = mpmath.sqrt(ai * gi)
        s = mpmath.ldexp(s, -1)
        q *= mpmath.power(g / a, s)
        if mpmath.fabs(a - g) <= mpmath.ldexp(g, -(mpmath.mp.prec - 1)):
            return q, a, k

    util.error("AGM iteration did not converge")

def sum_T(c, u):
    b0 = mpmath.mpf(0)
    b1 = mpmath.mpf(0)
    for ci in reversed(c):
        temp = b0
        b0 = 2 * u * b0 - b1 + ci
        b1 = temp

    return b0 - b1 * u

def sum_V(c, u):
    b0 = mpmath.mpf(0)
    b1 = mpmath.mpf(0)
    for ci in reversed(c):
        temp = b0
        b0 = 2 * u * b0 - b1 + ci
        b1 = temp

    return b0 - b1

class phi:
    def set_T_coef(self, eps):
        C = mpmath.fdiv(1, mpmath.sqrt(self.r * self.K2_pi) * self.mlogq)
        for N in range(1, 10000):
            e = mpmath.fmul(N, N)
            if mpmath.fdiv(C * mpmath.power(self.q, e), N) <= eps:
                self.c0 = (N + 1) * [mpmath.mpf(1)]
                for n in range(1, N + 1):
                    self.c0[n] = 2 * mpmath.power(self.q, mpmath.fmul(n, n))
                return

        util.error("iteration did not converge")

    def set_V_coef(self, eps):
        C = mpmath.mpf(2)
        C /= self.k * mpmath.power(self.r, mpmath.ldexp(1, -2))
        C /= mpmath.power(self.K2_pi, mpmath.ldexp(3, -1)) * self.mlogq
        for N in range(1, 10000):
            A = mpmath.fdiv(2 * N + 3, 2 * N + 1)
            Nh = mpmath.fadd(N, mpmath.ldexp(1, -1))
            B = mpmath.power(self.q2, Nh * Nh)
            if C * A * B <= eps:
                self.c1 = (N + 1) * [mpmath.mpf(1)]
                for n in range(1, N + 1):
                    e =  mpmath.fadd(mpmath.fmul(n, n), n)
                    self.c1[n] = mpmath.power(self.q2, e)
                return

        util.error("iteration did not converge")

    def __init__(self, r):
        self.r = mpmath.mpf(r)
        self.q, self.pi_2K, self.k = q_pi_2K_k(self.r)
        self.sr = mpmath.sqrt(self.r)
        self.q2 = self.q * self.q
        self.K2_pi = 1 / self.pi_2K
        self.rhoA = 1 / self.q
        self.mlogq = - mpmath.log(self.q)
        self.C1 = self.k * mpmath.power(self.r, mpmath.ldexp(3, -2))
        self.C1 *= mpmath.sqrt(self.q)
        self.C1 *= mpmath.power(self.K2_pi, mpmath.ldexp(3, -1))

        eps = mpmath.ldexp(1, - mpmath.mp.prec)
        self.set_T_coef(eps)
        self.set_V_coef(eps)

    @staticmethod
    def name():
        return "Phi"

    @staticmethod
    def latexname(r):
        return r"\Phi_{{{}}}".format(r)

    def value(self, u):
        return self.sr * sum_T(self.c0, u) / sum_T(self.c0, -u)

    def deriv1(self, u):
        T = sum_T(self.c0, -u)
        V = sum_V(self.c1, 1 - 2 * u * u)
        return self.C1 * V / (T * T)

    def hatrho(self):
        return self.rhoA

