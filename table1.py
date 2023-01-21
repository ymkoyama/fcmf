import util
import importlib.util
import mpmath
import mpmath_util

mpmath.mp.prec = 100
phi_list = ["Phi", "exp", "P2", "P1"]
phi_module_list = [importlib.import_module(phi) for phi in phi_list]
log2r_list = [-i for i in range(1, 21)]
def mpstr(x):
    return mpmath.nstr(x, n = 6, strip_zeros = False)

util.print_sys_version()
mpmath_util.print_version()
mpmath_util.print_prec()
print("phi", ",".join([phi_module.__name__ for phi_module in phi_module_list]))

title = "r"
for phi_module in phi_module_list:
    hatrho = "hatrho" + phi_module.phi.name()
    title += " " + hatrho + " " + hatrho + "^2"
print(title)
for log2r in log2r_list:
    r = mpmath.ldexp(1, log2r)
    hatrho_hatrho2_list = []
    for phi_module in phi_module_list:
        phi = phi_module.phi(r)
        hatrho = phi.hatrho()
        hatrho_hatrho2_list.append(hatrho)
        hatrho_hatrho2_list.append(hatrho * hatrho)

    print("$2^{" + str(log2r) + "}$", "&",
        " & ".join([mpstr(x) for x in hatrho_hatrho2_list]), "\\\\")

