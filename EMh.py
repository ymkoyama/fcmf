import argparse
import util
import importlib.util
import mpmath
import mpmath_util
import cmf
import time

total_time0 = time.time()

parser = argparse.ArgumentParser(
    description = "Save E_{M,h} and E_{M,h}/f(0)",
    formatter_class = argparse.RawTextHelpFormatter)

parser.add_argument("--cmf", default = "inverse_power",
    help = "Module name which defines CMF class")
parser.add_argument("--param", default = "1",
    help = "Parameter for CMF class")
parser.add_argument("--a", default = "0.5",
    help = "a for int_a^b exp(-xt) w(t) dt")
parser.add_argument("--b", default = "1",
    help = "See --a")
parser.add_argument("--method", default = "qf",
    help = "--method=qf computes with the quadratic form representation\n"
           "--method=det computes with the determinant representation")
parser.add_argument("--Mmin", type = int, default = 1,
    help = "Compute for Mmin <= M <= Mmax")
parser.add_argument("--Mmax", type = int, default = 4,
    help = "See --Mmin")
parser.add_argument("--N", type = int, default = 201,
    help = "Evaluate at h = e^linspace(log(hmin), log(hmax), N)\n"
           "if bound_i >= ymin or bound_{i-1} >= ymin or bound_{i+1} >= ymin")
parser.add_argument("--hmin", default = "0.1",
    help = "See --N")
parser.add_argument("--hmax", default = "10",
    help = "See --N")
parser.add_argument("--ymin", default = "1e-50",
    help = "See --N")
parser.add_argument("--prec", type = int, default = 200,
    help = "Set mpmath.mp.prec (bit)")

args = parser.parse_args()

if importlib.util.find_spec(args.cmf) is None:
    util.error(args.cmf + " can not be found")

if args.method != "qf" and args.method != "det":
    util.error("method=qf,det")

if args.Mmin <= 0 or args.Mmin > args.Mmax:
    util.error("1 <= Mmin <= Mmax")

if args.N <= 0:
    util.error("N >= 1")

mpmath.mp.prec = args.prec
cmf_module = importlib.import_module(args.cmf)
param = args.param
method = args.method
a = mpmath.mpf(args.a)
b = mpmath.mpf(args.b)
Mmin = args.Mmin
Mmax = args.Mmax
N = args.N
hmin = mpmath.mpf(args.hmin)
hmax = mpmath.mpf(args.hmax)
ymin = mpmath.mpf(args.ymin)
filename_zero_padding = 4

if a <= 0 or a >= b:
    util.error("0 < a < b")

if hmin <= 0 or hmin > hmax:
    util.error("0 < hmin <= hmax")

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
print("N", N)
print("hmin", hmin)
print("hmax", hmax)
print("ymin", ymin)

f = cmf_module.CMF(a, b, param)
f0 = f.value(mpmath.mpf(0))
print("f0", f0)
hr, F, A, B = cmf.hr_factors(f.a / f.b)
util.save_vector([f0], "f0.dat")
util.save_vector([f.a], "a.dat")
util.save_vector([f.b], "b.dat")
util.save_vector([hr], "hab.dat")
util.save_vector([hr / f.b], "hab_b.dat")
util.save_vector([F], "Fab.dat")
util.save_vector([A], "Aab.dat")
util.save_vector([B], "Bab.dat")

h = [mpmath.exp(lh) for lh in
    mpmath.linspace(mpmath.log(hmin), mpmath.log(hmax), N)]

for M in range(Mmin, Mmax + 1):
    print("M", M)
    hc = []
    EMh_f0_bound = []
    EMh_f0 = []
    EMh = []
    bound = [cmf.EMh_f0_bound(f.a, f.b, M, hi) for hi in h]
    for i in range(len(h)):
        flag = False
        if bound[i] >= ymin:
            flag = True
        elif i != len(h) - 1 and bound[i + 1] >= ymin:
            flag = True
        elif i != 0 and bound[i - 1] >= ymin:
            flag = True

        if flag:
            hc.append(h[i])
            EMh_f0_bound.append(bound[i])
            EMhi = cmf.EMh_determinant(f, M, h[i]) if method == "det" \
                else cmf.EMh_quadratic_form(f, M, h[i])
            EMh.append(EMhi)
            EMh_f0.append(EMhi / f0)

    suffix = "M" + str(M).zfill(filename_zero_padding) + ".dat"
    util.save_vector2(hc, EMh_f0_bound, "EMh_f0_bound_" + suffix)
    util.save_vector2(hc, EMh_f0, "EMh_f0_" + suffix)
    util.save_vector2(hc, EMh, "EMh_" + suffix)

dt = time.time() - total_time0
print("time_total(sec)", dt)

