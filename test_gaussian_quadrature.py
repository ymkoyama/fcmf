import mpmath
import gaussian_quadrature as gq

def gauss_legendre_quadrature_exact(M):
    u = []
    c = []

    #https://en.wikipedia.org/wiki/Gaussian_quadrature
    if M == 1:
        u = [mpmath.mpf(0)]
        c = [mpmath.mpf(2)]
    elif M == 2:
        a = mpmath.fdiv(1, mpmath.sqrt(3))
        u = [-a, a]
        c = 2 * [mpmath.mpf(1)]
    elif M == 3:
        a = mpmath.sqrt(mpmath.fdiv(3, 5))
        u = [-a, mpmath.mpf(0), a]
        c = [mpmath.fdiv(5, 9), mpmath.fdiv(8, 9), mpmath.fdiv(5, 9)]
    elif M == 4:
        a1 = mpmath.fdiv(3, 7)
        a2 = mpmath.fdiv(2, 7) * mpmath.sqrt(mpmath.fdiv(6, 5))
        u = [-mpmath.sqrt(a1 + a2), -mpmath.sqrt(a1 - a2),
            mpmath.sqrt(a1 - a2), mpmath.sqrt(a1 + a2)]
        c1 = mpmath.fdiv(18 - mpmath.sqrt(30), 36)
        c2 = mpmath.fdiv(18 + mpmath.sqrt(30), 36)
        c = [c1, c2, c2, c1]

    return u, c

def test_discretized_stieltjes():
    print("")
    mpmath.mp.prec = 100
    M = 4
    MDS = M
    u, c = gauss_legendre_quadrature_exact(MDS)
    alpha, beta = gq.discretized_stieltjes(M, u, c)
    alpha0, beta0 = gq.gauss_legendre_alpha_beta(M)

    print("alpha alpha0 alpha-alpha0")
    for i in range(len(alpha)):
        print(alpha[i], alpha0[i], alpha[i] - alpha0[i])
        assert mpmath.almosteq(alpha[i], alpha0[i])

    print("beta beta0 beta-beta0")
    for i in range(len(beta)):
        print(beta[i], beta0[i], beta[i] - beta0[i])
        assert mpmath.almosteq(beta[i], beta0[i])

def test_gauss_legendre_quadrature_golub_welsch():
    print("")
    mpmath.mp.prec = 100

    for M in range(1, 5):
        print("M", M)
        u0, c0 = gauss_legendre_quadrature_exact(M)
        assert len(u0) == M
        assert len(c0) == M

        u, c = gq.gauss_legendre_quadrature_golub_welsch(M)
        assert len(u) == M
        assert len(c) == M

        print("u u0 u-u0")
        for i in range(M):
            print(u[i], u0[i], u[i] - u0[i])
            assert mpmath.almosteq(u[i], u0[i])

        print("c c0 c-c0")
        for i in range(M):
            print(c[i], c0[i], c[i] - c0[i])
            assert mpmath.almosteq(c[i], c0[i])

def test_gauss_legendre_quadrature_mpmath():
    print("")
    mpmath.mp.prec = 100

    for exponent in range(3):
        M = 3 * 2 ** exponent
        um, cm = gq.gauss_legendre_quadrature_mpmath(M)
        assert len(um) == M
        assert len(cm) == M

        u, c = gq.gauss_legendre_quadrature_golub_welsch(M)
        assert len(u) == M
        assert len(c) == M

        print("exponent", exponent)
        print("M", M)

        print("u um u-um")
        for i in range(len(u)):
            print(u[i], um[i], u[i] - um[i])
            assert mpmath.almosteq(u[i], um[i])

        print("c cm c-cm")
        for i in range(len(c)):
            print(c[i], cm[i], c[i] - cm[i])
            assert mpmath.almosteq(c[i], cm[i])

