import mpmath
import cmf
import gaussian_quadrature as gq
import inverse_power
import P1

def test_inverse_power():
    mpmath.mp.prec = 100
    param_list = [mpmath.mpf("0.5"), mpmath.mpf("1"), mpmath.mpf("2")]
    a_list = [mpmath.ldexp(1, -1), mpmath.ldexp(1, -10), mpmath.mpf(2)]
    b_list = [mpmath.mpf(1),       mpmath.mpf(1),        mpmath.mpf(5)]
    n_list = range(5)
    x_list = mpmath.linspace(0, 10, 11)
    eps = mpmath.ldexp(1, -(mpmath.mp.prec - 10))
    print("\nmpmath.mp.prec", mpmath.mp.prec, "eps", eps)
    prec0 = 2 * mpmath.mp.prec

    def derivn_quad(f, x, n):
        g = lambda t: mpmath.power(-t, n) * mpmath.exp(- x * t) * f.w(t)
        return mpmath.quad(g, [f.a, f.b])

    for param in param_list:
        for a, b in zip(a_list, b_list):
            f = inverse_power.CMF(a, b, param)
            for n in n_list:
                print("param", param, "a", a, "b", b, "derivative", n)
                print("x derivn derivn_quad rel_error")
                for x in x_list:
                    v = f.derivn(x, n)
                    temp = mpmath.mp.prec
                    mpmath.mp.prec = prec0
                    vq = derivn_quad(f, x, n)
                    mpmath.mp.prec = temp
                    re = (v - vq) / vq
                    print(x, v, vq, mpmath.nstr(re))
                    assert mpmath.fabs(re) <= eps

def test_solve_y_newton():
    print("")
    #f(x) = pi*exp(-sqrt(2)*x)
    #f(x) = y -> x = log(pi/y)/sqrt(2)
    prec_list = [50, 100, 200, 500, 1000]
    print_log = True

    for prec in prec_list:
        mpmath.mp.prec = prec
        print("mpmath.mp.prec", mpmath.mp.prec)
        eps = mpmath.ldexp(1, -(mpmath.mp.prec - 10))
        y_list = mpmath.linspace(1, 3, 3)
        t = [mpmath.sqrt(2)]
        c = [mpmath.pi]
        f = cmf.ExponentialSum(t, c)
        for y in y_list:
            x0 = mpmath.log(mpmath.pi / y) / mpmath.sqrt(2)
            x = cmf.solve_y_newton(f, y, mpmath.mpf(0), print_log)
            re = (x - x0) / x
            print("y", y, "x", x, "x0", x0, "re", re)
            assert mpmath.fabs(re) <= eps

def test_solve_y_bisection_newton():
    print("")
    prec_list = [50, 100, 200, 500, 1000]
    print_log = True

    for prec in prec_list:
        mpmath.mp.prec = prec
        eps = mpmath.ldexp(1, -(mpmath.mp.prec - 10))
        param_list = [mpmath.mpf("0.5"), mpmath.mpf("1"), mpmath.mpf("2")]
        a_list = [mpmath.ldexp(1, -i) for i in range(1, 21)]
        b = mpmath.mpf(1)

        for param in param_list:
            for a in a_list:
                f = inverse_power.CMF(a, b, param)
                f0 = f.value(mpmath.mpf(0))
                fxmax0 = mpmath.ldexp(f0, -mpmath.mp.prec)
                xmax = cmf.solve_y_bisection_newton(f, fxmax0, print_log)
                print("mpmath.mp.prec", mpmath.mp.prec,
                    "param", param, "a", a, "b", b)
                fxmax = f.value(xmax)
                re = (fxmax - fxmax0) / fxmax0
                print("f(xmax)", fxmax, "f(0)*2^-prec", fxmax0, "re", mpmath.nstr(re))
                assert mpmath.fabs(re) <= eps

def test_hr_factors():
    print("")
    #(1-r)*exp(-hr) + exp(-(1-r)*hr) - r = 0
    #r=1/2 -> (1/2)*exp(-hr) + exp(-(1/2)*hr) - 1/2 = 0
    #      -> exp(-hr) + 2*exp(-(1/2)*hr) - 1 = 0
    #      -> (exp(-(1/2)*hr) + 1)^2 - 2 = 0
    #      -> exp(-(1/2)*hr) = sqrt(2) - 1
    #      -> exp(hr/2) = sqrt(2) + 1
    #      -> hr = 2*log(sqrt(2)+1)
    #      -> F = r*exp((1-r)*hr) = exp(hr/2)/2 = (sqrt(2) + 1)/2
    prec_list = [50, 100, 200, 500, 1000]
    print_log = True
    for prec in prec_list:
        mpmath.mp.prec = prec
        eps = mpmath.ldexp(1, -(mpmath.mp.prec - 1))

        hr0 = mpmath.ldexp(mpmath.log(mpmath.sqrt(2) + 1), 1)
        F0 = mpmath.ldexp(mpmath.sqrt(2) + 1, -1)
        r = mpmath.ldexp(1, -1)
        hr, F, A, B = cmf.hr_factors(r, print_log)

        re = (hr - hr0) / hr0
        print("prec", mpmath.mp.prec, "r", r, "hr", hr, "hr0", hr0, "re", re)
        assert mpmath.fabs(re) <= eps

        re = (F - F0) / F0
        print("prec", mpmath.mp.prec, "r", r, "F", F, "F0", F0, "re", re)
        assert mpmath.fabs(re) <= eps

        r_list = [mpmath.mpf("0.99"), mpmath.mpf("0.9")]
        r_list += [mpmath.ldexp(1, -n) for n in range(1, 21)]
        r_list += [mpmath.power(10, -10)]
        r_list += [mpmath.power(10, -100)]
        for r in r_list:
            hr, F, A, B = cmf.hr_factors(r, print_log)
            onemr = 1 - r
            G = onemr * mpmath.exp(-hr) + mpmath.exp(-onemr * hr) - r
            print("prec", mpmath.mp.prec, "r", r, "hr", hr, "G(hr)", G)
            assert mpmath.fabs(G) <= eps

def test_best_dF_dtcx():
    print("")
    mpmath.mp.prec = 100
    param = mpmath.mpf(1)
    a = mpmath.mpf(2)
    b = mpmath.mpf(5)
    f = inverse_power.CMF(a, b, param)

    for M in range(1, 4):
        t0 = mpmath.linspace(1, M, M)
        c0 = mpmath.linspace(M + 1, 2 * M, M)
        x0 = mpmath.linspace(2 * M + 1, 4 * M, 2 * M)
        dF_dtcx = cmf.best_dF_dtcx(f, t0, c0, x0)
        print("M", M)
        print("t", t0)
        print("c", c0)
        print("x", x0)
        print("dF_dtcx.rows", dF_dtcx.rows)
        print("dF_dtcx.cols", dF_dtcx.cols)
        assert dF_dtcx.rows == 4 * M
        assert dF_dtcx.cols == 4 * M

        print("i j Fij Fijdiff Fij-Fijdiff")
        for i in range(4 * M):
            if M == 1:
                Fi = lambda t, c, x1, x2: cmf.best_F(f, [t], [c], [x1, x2])[i]
            elif M == 2:
                Fi = lambda t1, t2, c1, c2, x1, x2, x3, x4: cmf.best_F(f, [t1, t2], [c1, c2], [x1, x2, x3, x4])[i]
            elif M == 3:
                Fi = lambda t1, t2, t3, c1, c2, c3, x1, x2, x3, x4, x5, x6: cmf.best_F(f, [t1, t2, t3], [c1, c2, c3], [x1, x2, x3, x4, x5, x6])[i]

            for j in range(4 * M):
                n = 4 * M * [0]
                n[j] = 1
                Fij = dF_dtcx[i, j]
                Fijdiff = mpmath.diff(Fi, t0 + c0 + x0, n)
                print(i, j, Fij, Fijdiff, mpmath.nstr(Fij - Fijdiff))
                assert mpmath.almosteq(Fij, Fijdiff)

def exponential_sum_chebyshev_polynomial(M):
    E0 = mpmath.mpf(1)

    if M == 1:
        #T5(u) = 5 u - 20 u^3 + 16 u^5
        #E(x) = T5(exp(-x)) = 5 exp(-x) + 16 exp(-5x) - 20 exp(-3x)
        #cos(y) = exp(-x) -> x = log(1/cos(y)), 0 <= y < pi/2
        #E(x) = cos(5y) -> E'(x) = 0 -> sin(5y) = 0
        #-> x = log(1/cos(n*pi/5)), 1 <= n <= 2
        f = cmf.ExponentialSum(
            [mpmath.mpf(1), mpmath.mpf(5)],
            [mpmath.mpf(5), mpmath.mpf(16)])
        t = [mpmath.mpf(3)]
        c = [mpmath.mpf(20)]
        x = [mpmath.log(1 / mpmath.cospi(mpmath.fdiv(n, 5))) for n in mpmath.linspace(1, 2, 2)]

        return f, t, c, x, E0

    elif M == 2:
        #T9(u) = 9 u - 120 u^3 + 432 u^5 - 576 u^7 + 256 u^9
        #E(x) = T9(exp(-x)) = 9 exp(-x) + 432 exp(-5x)  + 256 exp(-9x)
        # - 120 exp(-3x) - 576 exp(-7x)
        #cos(y) = exp(-x) -> x = log(1/cos(y)), 0 <= y < pi/2
        #E(x) = cos(9y) -> E'(x) = 0 -> sin(9y) = 0
        #-> x = log(1/cos(n*pi/9)), 1 <= n <= 4
        f = cmf.ExponentialSum(
            [mpmath.mpf(1), mpmath.mpf(5), mpmath.mpf(9)],
            [mpmath.mpf(9), mpmath.mpf(432), mpmath.mpf(256)])
        t = [mpmath.mpf(3), mpmath.mpf(7)]
        c = [mpmath.mpf(120), mpmath.mpf(576)]
        x = [mpmath.log(1 / mpmath.cospi(mpmath.fdiv(n, 9))) for n in mpmath.linspace(1, 4, 4)]

        return f, t, c, x, E0

    elif M == 3:
        #T13(u) = 13 u - 364 u^3 + 2912 u^5 - 9984 u^7
        # + 16640 u^9 - 13312 u^11 + 4096 u^13
        #E(x) = T13(exp(-x)) =
        # 13 exp(-x) + 2912 exp(-5x) + 16640 exp(-9x) + 4096 exp(-13x)
        # - 364 exp(-3x) - 9984 exp(-7x) - 13312 exp(-11x)
        #cos(y) = exp(-x) -> x = log(1/cos(y)), 0 <= y < pi/2
        #E(x) = cos(13y) -> E'(x) = 0 -> sin(13y) = 0
        #-> x = log(1/cos(n*pi/13)), 1 <= n <= 6
        f = cmf.ExponentialSum(
            [mpmath.mpf(t) for t in [1, 5, 9, 13]],
            [mpmath.mpf(c) for c in [13, 2912, 16640, 4096]])
        t = [mpmath.mpf(t) for t in [3, 7, 11]]
        c = [mpmath.mpf(c) for c in [364, 9984, 13312]]
        x = [mpmath.log(1 / mpmath.cospi(mpmath.fdiv(n, 13))) for n in mpmath.linspace(1, 6, 6)]

        return f, t, c, x, E0

def test_best_newton():
    mpmath.mp.prec = 100
    eps = mpmath.power(10, -20)
    print("\nmpmath.mp.prec", mpmath.mp.prec, "eps", mpmath.nstr(eps))
    for M in range(1, 4):
        print("M", M)
        f, t0, c0, x0, E0 = exponential_sum_chebyshev_polynomial(M)

        s = mpmath.mpf("1.01")
        t = [s * ti for ti in t0]
        c = [s * ci for ci in c0]
        x = [s * xi for xi in x0]

        t, c, x = cmf.best_newton(f, t, c, x, max_iter = 10, print_log = 10)

        print("t", "t0", "(t-t0)/t0")
        for nu in range(M):
            re = (t[nu] - t0[nu]) / t0[nu]
            print(t[nu], t0[nu], mpmath.nstr(re))
            assert mpmath.fabs(re) <= eps

        print("c", "c0", "(c-c0)/c0")
        for nu in range(M):
            re = (c[nu] - c0[nu]) / c0[nu]
            print(c[nu], c0[nu], mpmath.nstr(re))
            assert mpmath.fabs(re) <= eps

        print("x", "x0", "(x-x0)/x0")
        for i in range(2 * M):
            re = (x[i] - x0[i]) / x0[i]
            print(x[i], x0[i], mpmath.nstr(re))
            assert mpmath.fabs(re) <= eps

def test_best_remez_chebyshev_polynomial():
    print("")
    mpmath.mp.prec = 100
    eps = mpmath.power(10, -20)
    print("mpmath.mp.prec", mpmath.mp.prec, "eps", mpmath.nstr(eps))
    remez_max_scale = mpmath.mpf("0.1")
    remez_armijo = mpmath.mpf("0.1")
    remez_eps = mpmath.power(10, -5)
    remez_iter = 100

    for M in range(1, 4):
        f, t0, c0, x0, E0 = exponential_sum_chebyshev_polynomial(M)
        a = f.t[0]
        b = f.t[M]
        hr = cmf.hr(a / b)
        h = hr / b
        y = [mpmath.exp(- h * t) for t in f.t]
        C = [c * (1 + mpmath.exp(- h * t)) for t, c in zip(f.t, f.c)]
        alpha, beta = gq.discretized_stieltjes(M, y, C)
        print("M", M)
        print("a", a)
        print("b", b)
        print("a/b", a / b)
        print("h", h)

        t, c, x = cmf.init_remez_gauss(alpha, beta, h)
        assert len(t) == M
        assert len(c) == M
        assert len(x) == 2 * M

        print("t c")
        for i in range(M):
            print(t[i], c[i])

        print("x i*h x-i*h")
        for i in range(len(x)):
            print(x[i], (i + 1) * h, x[i] - (i + 1) * h)
            assert mpmath.almosteq(x[i], (i + 1) * h)

        fM = cmf.ExponentialSum(t, c)
        Ex = [f.value(xi) - fM.value(xi) for xi in [mpmath.mpf(0)] + x]
        print("xi", "E(xi)", "(E(x{i-1})+E(xi))/E(xi)")
        for i in range(len(x)):
            print(x[i], Ex[i + 1], mpmath.nstr((Ex[i] + Ex[i + 1]) / Ex[i + 1]))
            assert mpmath.fabs((Ex[i] + Ex[i + 1]) / Ex[i + 1]) <= eps

        t, c, x, iteration = cmf.best_remez(f, t, c, x,
            remez_max_scale, remez_armijo, remez_eps, remez_iter, print_log = 10)

        print("iteration", iteration)

        print("t", "t0", "(t-t0)/t0")
        for i in range(M):
            re = (t[i] - t0[i]) / t0[i]
            print(t[i], t0[i], mpmath.nstr(re))
            assert mpmath.fabs(re) <= remez_eps

        print("c", "c0", "(c-c0)/c0")
        for i in range(M):
            re = (c[i] - c0[i]) / c0[i]
            print(c[i], c0[i], mpmath.nstr(re))
            assert mpmath.fabs(re) <= remez_eps

        print("x", "x0", "(x-x0)/x0")
        for i in range(2 * M):
            re = (x[i] - x0[i]) / x0[i]
            print(x[i], x0[i], mpmath.nstr(re))
            assert mpmath.fabs(re) <= remez_eps

def test_best_remez_inverse_power1():
    print("")
    mpmath.mp.prec = 100
    eps = mpmath.power(10, -20)
    remez_max_scale = mpmath.mpf("0.1")
    remez_armijo = mpmath.mpf("0.1")
    remez_eps = mpmath.power(10, -5)
    remez_iter = 100
    MDS = 192
    Mmax = 3
    param = mpmath.mpf(1)
    a = mpmath.mpf(2)
    b = mpmath.mpf(5)
    f = inverse_power.CMF(a, b, param)
    hr = cmf.hr(a / b)
    h = hr / b
    print("a", a)
    print("b", b)
    print("a/b", a / b)
    print("h", h)

    print(MDS, "point Gauss-Legendre quadrature")
    uGL, cGL = gq.gauss_legendre_quadrature_mpmath(MDS)
    init_remez_phi = P1.phi(a / b)
    yDS, dDS = cmf.init_remez_yDS_dDS_gauss_legendre(
        f, h, init_remez_phi, uGL, cGL)
    alpha, beta = gq.discretized_stieltjes(Mmax, yDS, dDS)

    for M in range(1, Mmax + 1):
        print("M", M)
        t, c, x = cmf.init_remez_gauss(alpha[0 : M], beta[0 : M], h)
        assert len(t) == M
        assert len(c) == M
        assert len(x) == 2 * M

        print("t c")
        for i in range(M):
            print(t[i], c[i])

        print("x i*h x-i*h")
        for i in range(len(x)):
            print(x[i], (i + 1) * h, x[i] - (i + 1) * h)
            assert mpmath.almosteq(x[i], (i + 1) * h)

        fM = cmf.ExponentialSum(t, c)
        E = [f.value(xi) - fM.value(xi) for xi in [mpmath.mpf(0)] + x]
        print("xi", "E(xi)", "(E(x{i-1})+E(xi))/E(xi)")
        for i in range(len(x)):
            print(x[i], E[i + 1], mpmath.nstr((E[i] + E[i + 1]) / E[i + 1]))
            assert mpmath.fabs((E[i] + E[i + 1]) / E[i + 1]) <= eps

        t, c, x, iteration = cmf.best_remez(f, t, c, x,
            remez_max_scale, remez_armijo, remez_eps, remez_iter, print_log = 10)

        print("iteration", iteration)

        fM = cmf.ExponentialSum(t, c)
        E = [f.value(xi) - fM.value(xi) for xi in [mpmath.mpf(0)] + x]
        print("xi", "E(xi)", "(E(x{i-1})+E(xi))/E(xi)")
        for i in range(len(x)):
            print(x[i], E[i + 1], mpmath.nstr((E[i] + E[i + 1]) / E[i + 1]))
            assert mpmath.fabs((E[i] + E[i + 1]) / E[i + 1]) <= remez_eps


        E = [f.value(xi) - fM.value(xi) for xi in x]
        dE = [f.deriv1(xi) - fM.deriv1(xi) for xi in x]
        print("xi", "E'(xi)", "xi*E'(xi)/E(xi)")
        for i in range(len(x)):
            print(x[i], dE[i], mpmath.nstr(x[i] * dE[i] / E[i]))
            assert mpmath.fabs(x[i] * dE[i] / E[i]) <= remez_eps

def test_fcmf_gauss():
    print("")
    mpmath.mp.prec = 100
    param = mpmath.mpf(1)
    a = mpmath.mpf(2)
    b = mpmath.mpf(5)
    f = inverse_power.CMF(a, b, param)
    phi = P1.phi(a / b)
    Mmax = 4
    MDS = 48
    uGL, cGL = gq.gauss_legendre_quadrature_mpmath(MDS)
    uDS, cDS = cmf.uDS_cDS_gauss_legendre(f, phi, uGL, cGL)
    alpha, beta = gq.discretized_stieltjes(Mmax, uDS, cDS)

    for M in range(1, Mmax + 1):
        print("M", M)
        u, c = gq.golub_welsch(alpha[0 : M], beta[0 : M])
        assert len(u) == M
        assert len(c) == M

        uGL, cGL = gq.gauss_legendre_quadrature_golub_welsch(M)
        assert len(uGL) == M
        assert len(cGL) == M
        u0 = uGL
        c0 = [cGLi * (b - a) / 2 for cGLi in cGL]

        print("u u0 u-u0")
        for i in range(M):
            print(u[i], u0[i], u[i] - u0[i])
            assert mpmath.almosteq(u[i], u0[i])

        print("c c0 c-c0")
        for i in range(M):
            print(c[i], c0[i], c[i] - c0[i])
            assert mpmath.almosteq(c[i], c0[i])

def test_EMh():
    print("")
    mpmath.mp.prec = 100
    param = mpmath.mpf(1)
    a = mpmath.mpf(1)
    b = mpmath.mpf(2)
    f = inverse_power.CMF(a, b, param)

    print("M", "h", "EMh_quad", "EMh_det", "EMh_quad-EMh_det")
    for M in range(1, 4):
        for h in mpmath.linspace(1, 3, 3):
            EMh_quad = cmf.EMh_quadratic_form(f, M, h)
            EMh_det = cmf.EMh_determinant(f, M, h)
            print(M, h, EMh_quad, EMh_det, EMh_quad - EMh_det)
            assert mpmath.almosteq(EMh_quad, EMh_det)

def test_find_max_abs_error():
    mpmath.mp.prec = 100
    cutoff = mpmath.mpf("0.01")
    E = lambda x: mpmath.sin(x) - x * mpmath.cos(x)
    dE = lambda x: x * mpmath.sin(x)

    x = mpmath.linspace(0, mpmath.pi / 2, 101)
    max_abs_E = cmf.find_max_abs_error(E, dE, x, cutoff)
    assert mpmath.almosteq(max_abs_E, mpmath.mpf(1))

    x = mpmath.linspace(0, 4, 101)
    max_abs_E = cmf.find_max_abs_error(E, dE, x, cutoff)
    assert mpmath.almosteq(max_abs_E, mpmath.pi)

    x = mpmath.linspace(0, 8, 101)
    max_abs_E = cmf.find_max_abs_error(E, dE, x, cutoff)
    assert mpmath.almosteq(max_abs_E, 2 * mpmath.pi)

    x = mpmath.linspace(0, 12, 101)
    max_abs_E = cmf.find_max_abs_error(E, dE, x, cutoff)
    assert mpmath.almosteq(max_abs_E, -E(mpmath.mpf(12)))

    x = mpmath.linspace(0, 13, 101)
    max_abs_E = cmf.find_max_abs_error(E, dE, x, cutoff)
    assert mpmath.almosteq(max_abs_E, 4 * mpmath.pi)

    E = lambda x: - mpmath.sin(x) + (x + mpmath.pi) * mpmath.cos(x)
    dE = lambda x: - (x + mpmath.pi) * mpmath.sin(x)

    x = mpmath.linspace(0, mpmath.pi / 2, 101)
    max_abs_E = cmf.find_max_abs_error(E, dE, x, cutoff)
    assert mpmath.almosteq(max_abs_E, mpmath.pi)

    x = mpmath.linspace(0, 4, 101)
    max_abs_E = cmf.find_max_abs_error(E, dE, x, cutoff)
    assert mpmath.almosteq(max_abs_E, 2 * mpmath.pi)

    x = mpmath.linspace(0, 8, 101)
    max_abs_E = cmf.find_max_abs_error(E, dE, x, cutoff)
    assert mpmath.almosteq(max_abs_E, 3 * mpmath.pi)

    x = mpmath.linspace(0, 9, 101)
    max_abs_E = cmf.find_max_abs_error(E, dE, x, cutoff)
    assert mpmath.almosteq(max_abs_E, - E(mpmath.mpf(9)))

    x = mpmath.linspace(0, 10, 101)
    max_abs_E = cmf.find_max_abs_error(E, dE, x, cutoff)
    assert mpmath.almosteq(max_abs_E, 4 * mpmath.pi)

def test_phi_basis():
    print("")
    mpmath.mp.prec = 100
    eps = mpmath.power(10, -10)
    x = [mpmath.mpf(0)] + [mpmath.power(2, i) for i in mpmath.linspace(-10, 10, 21)]
    nmax = 4

    for log2r in [-1, -10]:
        r = mpmath.ldexp(1, log2r)
        print("r", r)
        phi = P1.phi(r)
        M = cmf.phi_basis_M(phi, nmax, eps)
        hatrho = phi.hatrho()
        print("hatrho", hatrho)
        print("nmax", nmax)
        print("M", M)
        for n in range(0, nmax + 1):
            s = mpmath.power(hatrho, n)
            print("phi_basis", "n", n, "s", s)
            print("x F F0 F-F0 s*(F-F0)")
            for xi in x:
                F = cmf.phi_basis(phi, n, xi, M)
                F0 = phi.basis(n, xi)
                print(xi, F, F0, mpmath.nstr(F - F0), mpmath.nstr(s * (F - F0)))
                assert s * mpmath.fabs(F - F0) <= eps

