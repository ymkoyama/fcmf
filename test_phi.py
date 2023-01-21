import importlib.util
import mpmath
import Phi

def test_phi():
    print("")
    mpmath.mp.prec = 100
    phi_list = ["Phi", "exp", "P2", "P1", "R0_1"]
    eps = mpmath.ldexp(1, -(mpmath.mp.prec - 10))
    print("mpmath.mp.prec", mpmath.mp.prec, "eps", mpmath.nstr(eps))
    phi_module_list = [importlib.import_module(phi) for phi in phi_list]
    for phi_module in phi_module_list:
        for log2r in [-1, -10, -20]:
            r = mpmath.ldexp(1, log2r)
            phi = phi_module.phi(r)
            print("#phi.name", phi.name(), "r", r)

            #test phi(-1)=r
            one = mpmath.mpf(1)
            v = phi.value(-one)
            re = (v - r) / r
            print("phi(-1)", v, r, "re", mpmath.nstr(re))
            assert mpmath.fabs(re) <= eps

            #test phi(1)=1
            v = phi.value(one)
            re = (v - one) / one
            print("phi(1)", v, one, "re", mpmath.nstr(re))
            assert mpmath.fabs(re) <= eps

            #test strictly increasing
            u = mpmath.linspace(-1, 1, 21)
            v = [phi.value(ui) for ui in u]
            for i in range(len(u) - 1):
                assert v[i] < v[i + 1]

            #test phi.deriv1
            print("u phi.deriv1 mpmath.diff1 rel_error")
            for ui in u:
                deriv1 = phi.deriv1(ui)
                diff1 = mpmath.diff(lambda u: phi.value(u), ui, n = 1)
                re = (deriv1 - diff1) / diff1
                print(ui, deriv1, diff1, mpmath.nstr(re))
                assert mpmath.fabs(re) <= eps

def test_q_pi_2K_k():
    print("")
    prec0 = 2000
    r_list = [mpmath.mpf(s) for s in ["0.99", "0.9", "0.8", "0.7", "0.6"]]
    r_list += [mpmath.ldexp(1, -n) for n in range(1, 101)]

    for prec in [50, 100, 200, 500, 1000]:
        mpmath.mp.prec = prec
        eps = mpmath.ldexp(1, -(mpmath.mp.prec - 10))
        print("mpmath.mp.prec", mpmath.mp.prec, "eps", mpmath.nstr(eps))

        for r in r_list:
            q, pi_2K, k = Phi.q_pi_2K_k(r)

            mpmath.mp.prec = prec0
            m0 = 1 - r * r
            q0 = mpmath.qfrom(m = m0)
            K0 = mpmath.ellipk(m0)
            pi_2K0 = mpmath.pi / (2 * K0)
            mpmath.mp.prec = prec

            re = (q - q0) / q0
            print("r", r, "q", q, "q0", q0, "re", mpmath.nstr(re))
            assert mpmath.fabs(re) <= eps

            re = (pi_2K - pi_2K0) / pi_2K0
            print("r", r, "pi_2K", pi_2K, "pi_2K0", pi_2K0,
                "re", mpmath.nstr(re))
            assert mpmath.fabs(re) <= eps

def test_Phi():
    print("")
    prec0 = 2000
    for prec in [50, 100, 200, 500, 1000]:
        mpmath.mp.prec = prec
        eps = mpmath.ldexp(1, - (mpmath.mp.prec - 10))
        print("mpmath.mp.prec", mpmath.mp.prec, "eps", mpmath.nstr(eps))
        for log2r in [-1, -10, -20]:
            r = mpmath.ldexp(1, log2r)
            phi = Phi.phi(r)
            print("r", r, "q", phi.q)

            #analytical values
            mpmath.mp.prec = prec0
            m0 = 1 - r * r
            K0 = mpmath.ellipk(m0)
            K_pi0 = K0 / mpmath.pi
            uana = mpmath.linspace(-1, 1, 3)
            vana = [r, mpmath.sqrt(r), mpmath.mpf(1)]
            dana = [mpmath.mpf(0)] * 3
            dana[1] = mpmath.sqrt(r) * (1 - r) * K_pi0
            dana[2] = m0 * K_pi0 * K_pi0
            dana[0] = r * dana[2]
            mpmath.mp.prec = prec

            for u, v, d in zip(uana, vana, dana):
                value = phi.value(u)
                deriv1 = phi.deriv1(u)

                re = (value - v) / v
                print("u", u, "phi.value", value, "valueana", v,
                    "re", mpmath.nstr(re))
                assert mpmath.fabs(re) <= eps

                re = (deriv1 - d) / d
                print("u", u, "phi.deriv1", deriv1, "deriv1ana", d,
                    "re", mpmath.nstr(re))
                assert mpmath.fabs(re) <= eps

            for u in mpmath.linspace(mpmath.mpf("-0.9"), mpmath.mpf("0.9"), 19):
                mpmath.mp.prec = prec0
                k2 = 1 - r * r
                K = mpmath.ellipk(k2)
                K_pi = K / mpmath.pi
                value0 = mpmath.ellipfun("dn")(K_pi * mpmath.acos(u), k2)
                deriv10 = mpmath.diff(lambda x: mpmath.ellipfun("dn")(K_pi * mpmath.acos(x), k2), u)
                mpmath.mp.prec = prec

                value = phi.value(u)
                re = (value - value0) / value0
                print("u", u, "phi.value", value, "dnacos", value0,
                    "re", mpmath.nstr(re))
                assert mpmath.fabs(re) <= eps

                deriv1 = phi.deriv1(u)
                re = (deriv1 - deriv10) / deriv10
                print("u", u, "phi.deriv1", deriv1, "diff_dnacos", deriv10,
                    "re", mpmath.nstr(re))
                assert mpmath.fabs(re) <= eps

