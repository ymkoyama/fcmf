import argparse
import util
import mpmath_util
import os
import importlib.util
import mpmath
import gaussian_quadrature as gq
import cmf
import time

total_time0 = time.time()

parser = argparse.ArgumentParser(
    description =
        "Save parameters of exponential sum approximation\n"
        "of the finite completely monotone function\n"
        "f(x) = int_a^b exp(-x * t) w(t) dt, 0 < a < b < inf.\n"
        "f(x) = sum_{i=1}^M ci * exp(-x * ti) + E_{a,b,M}(x)\n",
    formatter_class = argparse.RawTextHelpFormatter)

method_candidate = ["best", "gauss", "best_init"]

parser.add_argument("--cmf", default = "inverse_power",
    help = "Module name which defines CMF class")
parser.add_argument("--param", default = "1",
    help = "Parameter for CMF class")
parser.add_argument("--a", default = "0.5",
    help = "a for int_a^b exp(-xt) w(t) dt")
parser.add_argument("--b", default = "1",
    help = "See --a")
parser.add_argument("--method", default = "best",
    help = "A method to compute exponential sum approximation.\n"
           "--method=best best exponential sum approximation\n"
           "--method=gauss Gaussian quadrature\n"
           "--method=best_init only performs initialization of best exponential sum")
parser.add_argument("--Mmax", type = int, default = 4,
    help = "Compute nodes and weights for Mmin <= M <= Mmax")
parser.add_argument("--Mmin", type = int, default = 1,
    help = "See --Mmax")
parser.add_argument("--phi", default = "Phi",
    help = "Module name which defines phi class for Gaussian quadrature\n"
           "to compute exponential sum approximation.\n"
           "Multiple phi can be selected by separating comma:\n"
           "e.g. --phi=Phi,exp,P2,P1,R0_1")
parser.add_argument("--nx_M", type = int, default = 500,
    help = "x is sampled for 0, e^linspace(log(xmin),log(xmax),nx_M*M+1)\n"
           "xmin = xmintmean / tmean\n"
           "tmean = -f'(0)/f(0)\n"
           "xmax is determined to satisfy f(xmax)/f(0) = 2^-mpmath.mp.prec")
parser.add_argument("--xmintmean", default = "1e-6",
    help = "See --nx_M")
parser.add_argument("--basis_nmax_2Mmax", type = int, default = 4,
    help = "nmax = 2 * Mmax * basis_nmax_2Mmax")
parser.add_argument("--remez_max_scale", default = "0.1",
    help = "For the update of positions in the Remez algorithm,\n"
           "the increment of the position is limited by\n"
           "remez_max_scale * (difference of two positions)")
parser.add_argument("--remez_armijo", default = "0.4",
    help = "For the update of positions in the Remez algorithm,\n"
           "dx is determined to satisfy Armijo's condition\n"
           "E(x) < 0: E(x + dx) - E(x) <= armijo * E'(x) * dx\n"
           "E(x) > 0: E(x + dx) - E(x) >= armijo * E'(x) * dx")
parser.add_argument("--remez_eps", default = "1e-10",
    help = "Iteration of Remez algorithm is repeated to satisfy\n"
           "max|(E(x{i-1})+E(xi))/E(xi)| <= remez_eps\n"
           "max|xi*E'(xi)/E(xi)| <= remez_eps")
parser.add_argument("--remez_iter", type = int, default = 100000,
    help = "The maximum number of iteration for Remez algorithm")
parser.add_argument("--best_newton", type = int, default = 5,
    help = "The number of iteration for Newton method after Remez algorithm")
parser.add_argument("--print_log", type = int, default = 1,
    help = "Print option of Remez algorithm\n"
           "0: no, 1: summary, 2: summary and parameters")
parser.add_argument("--MDS", type = int, default = 768,
    help = "The node number of discretized Stieltjes procedure\n"
           "for Gaussian quadrature. Current implementation requires"
           "MDS=3*2^m, i.e. 3,6,12,24,48,96,192,384,768,1536,3072,6144,...")
parser.add_argument("--prec", type = int, default = 200,
    help = "Set mpmath.mp.prec (bit)")

args = parser.parse_args()

if importlib.util.find_spec(args.cmf) is None:
    util.error(args.cmf + " can not be found")

if args.Mmin <= 0 or args.Mmin > args.Mmax:
    util.error("1 <= Mmin <= Mmax")

if args.method not in method_candidate:
    util.error("method=" + args.method + " is not found.\n"
       "method=" + ",".join(method_candidate))
method = args.method

phi_list = args.phi.split(",")
if len(phi_list) == 0:
    util.error("phi is empty")
for phi in phi_list:
    if importlib.util.find_spec(phi) is None:
        util.error("phi=" + phi + " can not be found")
phi_module_list = [importlib.import_module(phi) for phi in phi_list]

if args.remez_iter < 0:
    util.error("remez_iter >= 0")

if args.best_newton < 0:
    util.error("best_newton >= 0")

if args.print_log < 0:
    util.error("print_log >= 0")

if args.MDS <= 0:
    util.error("MDS >= 1")

def mpmath_eval(s):
    return mpmath.mpf(s) if s.isnumeric() else eval(s)

mpmath.mp.prec = args.prec
cmf_module = importlib.import_module(args.cmf)
param = mpmath_eval(args.param)
a = mpmath_eval(args.a)
b = mpmath_eval(args.b)
Mmin = args.Mmin
Mmax = args.Mmax
xmintmean = mpmath.mpf(args.xmintmean)
nx_M = args.nx_M
basis_nmax_2Mmax = args.basis_nmax_2Mmax
remez_max_scale = mpmath.mpf(args.remez_max_scale)
remez_armijo = mpmath.mpf(args.remez_armijo)
remez_iter = args.remez_iter
remez_eps = mpmath.mpf(args.remez_eps)
best_newton = args.best_newton
print_log = args.print_log
MDS = args.MDS
filename_zero_padding = 4
max_abs_E_cutoff = mpmath.mpf("0.01")

if a <= 0 or a >= b:
    util.error("0 < a < b")

if remez_max_scale <= 0:
    util.error("remez_max_scale > 0")

if remez_armijo <= 0:
    util.error("remez_armijo > 0")

if remez_eps <= 0:
    util.error("remez_eps > 0")

if xmintmean <= 0:
    util.error("xmintmean > 0")

if nx_M <= 0:
    util.error("nx_M >= 1")

if basis_nmax_2Mmax <= 0:
    util.error("basis_nmax_2Mmax >= 1")

def save_error(x, f, t, c, filename):
    fM = cmf.ExponentialSum(t, c)
    out = open(filename, "w")
    for i in range(len(x)):
        fi = f.value(x[i])
        fMi = fM.value(x[i])
        Ei = fi - fMi
        print(x[i], fi, fMi, Ei, file = out)
    out.close()
    print(filename)

def save_max_abs_error(x, f, t, c, filename):
    fM = cmf.ExponentialSum(t, c)
    E = lambda x: f.value(x) - fM.value(x)
    dE = lambda x: f.deriv1(x) - fM.deriv1(x)
    max_abs_E = cmf.find_max_abs_error(E, dE, x, max_abs_E_cutoff)
    util.save_vector([max_abs_E], filename)
    print("max_abs_E", mpmath.nstr(max_abs_E))
    return max_abs_E

util.print_sys_version()
mpmath_util.print_version()
mpmath_util.print_prec()
print("cmf", cmf_module.__name__)
print("param", param)
print("a", a)
print("b", b)
print("method", method)
print("Mmin", Mmin)
print("Mmax", Mmax)
print("phi", ",".join([phi_module.__name__ for phi_module in phi_module_list]))
print("xmintmean", xmintmean)
print("nx_M", nx_M)
print("basis_nmax_2Mmax", basis_nmax_2Mmax)
print("remez_max_scale", remez_max_scale)
print("remez_armijo", remez_armijo)
print("remez_iter", remez_iter)
print("remez_eps", remez_eps)
print("best_newton", best_newton)
print("print_log", print_log)
print("MDS", MDS)

f = cmf_module.CMF(a, b, param)
f0 = f.value(mpmath.mpf(0))
util.save_vector([f0], "f0.dat")
util.save_vector([f.a], "a.dat")
util.save_vector([f.b], "b.dat")

print("f(0)", mpmath.nstr(f0))
tmean = cmf.tmean(f)
print("tmean", mpmath.nstr(tmean))
print("1/tmean", mpmath.nstr(1 / tmean))
xmin = xmintmean / tmean
print("xmin", mpmath.nstr(xmin))
fxmax = mpmath.ldexp(f0, -mpmath.mp.prec)
xmax = cmf.solve_y_bisection_newton(f, fxmax)
print("xmax", mpmath.nstr(xmax))
print("f(xmax)/f(0)", mpmath.nstr(f.value(xmax) / f0))

time0 = time.time()
print(MDS, "point Gauss-Legendre quadrature", end = "", flush = True)
uGL, cGL = gq.gauss_legendre_quadrature_mpmath(MDS)
print(" time(sec)", time.time() - time0)

if method == "best" or method == "best_init":
    print("#method=", method, sep = "")

    hr = cmf.hr(a / b)
    h = hr / b
    print("hr", mpmath.nstr(hr))
    print("h", mpmath.nstr(h))
    util.save_vector([h], "best_init_h.dat")

    for phi_module in phi_module_list:
        phi = phi_module.phi(a / b)
        print("init_remez_alpha_beta", end = "", flush = True)
        time0 = time.time()
        yDS, dDS = cmf.init_remez_yDS_dDS_gauss_legendre(f, h, phi, uGL, cGL)
        alpha, beta = gq.discretized_stieltjes(Mmax, yDS, dDS)
        print(" time(sec)", time.time() - time0)
        prefix_init = "best_init_" + phi.name()
        suffix = "Mmax" + str(Mmax).zfill(filename_zero_padding) + ".dat"
        util.save_vector(alpha, prefix_init + "_alpha_" + suffix)
        util.save_vector(beta, prefix_init + "_beta_" + suffix)

        for M in range(Mmin, Mmax + 1):
            print("M", M, "init_remez_gauss", end = "", flush = True)
            time0 = time.time()
            t0, c0, x0 = cmf.init_remez_gauss(alpha[0 : M], beta[0 : M], h)
            E0 = f.value(mpmath.mpf(0)) - mpmath.fsum(c0)
            print(" time(sec)", time.time() - time0, "E0", mpmath.nstr(E0))
            suffix = "M" + str(M).zfill(filename_zero_padding) + ".dat"
            util.save_vector([E0], prefix_init + "_E0_" + suffix)
            util.save_vector(t0, prefix_init + "_t_" + suffix)
            util.save_vector(c0, prefix_init + "_c_" + suffix)
            util.save_vector(x0, prefix_init + "_x_" + suffix)

            save_x = mpmath.linspace(0, 2 * M * h, M * nx_M, endpoint = False)
            if 2 * M * h < xmax:
                save_x += [mpmath.exp(lx) for lx in
                    mpmath.linspace(mpmath.log(2 * M * h), mpmath.log(xmax),
                        nx_M + 1)]
            else:
                save_x += [2 * M * h]
            save_error(save_x, f, t0, c0, prefix_init + "_E_" + suffix)


            if method == "best_init":
                continue

            print("M", M, "best_remez")
            time0 = time.time()
            t, c, xst, iteration = cmf.best_remez(f, t0, c0, x0,
                remez_max_scale, remez_armijo, remez_eps, remez_iter, print_log)
            print("M", M, "best_remez", "time(sec)", time.time() - time0,
                "iteration", iteration)

            if best_newton >= 1:
                print("M", M, "best_newton")
                time0 = time.time()
                t, c, xst = cmf.best_newton(f, t, c, xst, best_newton, print_log)
                print("M", M, "best_newton",
                    "time(sec)", time.time() - time0, "iteration", best_newton)

            prefix = "best_" + phi.name()
            util.save_vector(t, prefix + "_t_" + suffix)
            util.save_vector(c, prefix + "_c_" + suffix)
            util.save_vector(xst, prefix + "_xst_" + suffix)
            x = [mpmath.mpf(0)] + [mpmath.exp(lx) for lx
                in mpmath.linspace(mpmath.log(xmin), mpmath.log(xmax),
                    nx_M * M + 1)]
            save_error(x, f, t, c, prefix + "_E_" + suffix)
            save_max_abs_error(x, f, t, c, prefix + "_max_abs_E_" + suffix)

        if method == "best":
            break

elif method == "gauss":
    for phi_module in phi_module_list:
        phi = phi_module.phi(a / b)
        prefix = method + "_" + phi.name()

        print("#method", method, "phi", phi.name())

        print("gauss_alpha_beta", end = "", flush = True)
        time0 = time.time()
        uDS, cDS = cmf.uDS_cDS_gauss_legendre(f, phi, uGL, cGL)
        alpha, beta = gq.discretized_stieltjes(Mmax, uDS, cDS)
        print(" time(sec)", time.time() - time0)
        suffix = "Mmax" + str(Mmax).zfill(filename_zero_padding) + ".dat"
        util.save_vector(alpha, prefix + "_alpha_" + suffix)
        util.save_vector(beta, prefix + "_beta_" + suffix)

        nmax = 2 * Mmax * basis_nmax_2Mmax
        n = list(range(0, nmax + 1))
        basis_sigma = [cmf.phi_basis_sigma(ni, uDS, cDS) for ni in n]
        nmaxstr = "nmax" + str(nmax).zfill(filename_zero_padding)
        suffix = nmaxstr + ".dat"
        util.save_vector(basis_sigma, prefix + "_basis_sigma_" + suffix)

        for M in range(Mmin, Mmax + 1):
            print("M", M, "golub_welsch", end = "", flush = True)
            time0 = time.time()
            u, c = gq.golub_welsch(alpha[0 : M], beta[0 : M])
            print(" time(sec)", time.time() - time0)

            Mstr = "M" + str(M).zfill(filename_zero_padding)
            suffix = Mstr + ".dat"
            util.save_vector(u, prefix + "_u_" + suffix)
            t = [f.b * phi.value(ui) for ui in u]
            util.save_vector(t, prefix + "_t_" + suffix)
            util.save_vector(c, prefix + "_c_" + suffix)

            x = [mpmath.mpf(0)] + [mpmath.exp(lx) for lx
                in mpmath.linspace(mpmath.log(xmin), mpmath.log(xmax),
                    nx_M * M + 1)]
            save_error(x, f, t, c, prefix + "_E_" + suffix)
            save_max_abs_error(x, f, t, c, prefix + "_max_abs_E_" + suffix)

            basis_sigma_gauss = [cmf.phi_basis_sigma(ni, u, c) for ni in n]
            basis_eps = [bs - bsg for bs, bsg in
                zip(basis_sigma, basis_sigma_gauss)]
            suffix = Mstr + "_" + nmaxstr + ".dat"
            util.save_vector(basis_eps, prefix + "_basis_eps_" + suffix)

print("time_total(sec)", time.time() - total_time0)

