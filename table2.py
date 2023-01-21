import util
import mpmath
import mpmath_util
import cmf

mpmath.mp.prec = 100
log2r_list = [-i for i in range(1, 21)]
def mpstr(x):
    return mpmath.nstr(x, n = 8, strip_zeros = False)

util.print_sys_version()
mpmath_util.print_version()
mpmath_util.print_prec()

print("r", "hr", "F=r*exp((1-r)*hr)", "A=(sqrt(F)+1)^2",
    "B=((sqrt(F)+1)/(sqrt(F)-1))^2")
for log2r in log2r_list:
    r = mpmath.ldexp(1, log2r)
    hr, F, A, B = cmf.hr_factors(r)
    print("$2^{" + str(log2r) + "}$", "&",
        " & ".join([mpstr(x) for x in [hr, F, A, B]]), "\\\\")

