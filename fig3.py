import util
import mpmath
import mpmath_util
import matplotlib.pyplot as plt
import matplotlib
import Phi
import exp
import P2
import P1

mpmath.mp.prec = 100
N = 201

util.print_sys_version()
mpmath_util.print_version()
mpmath_util.print_prec()
print("N", N)

def plot_hatrho2(ax, r):
    for pm, i in zip(phi_module_list, range(len(phi_module_list))):
        hatrho2 = [mpmath.power(pm.phi(ri).hatrho(), 2) for ri in r]
        name = pm.phi.latexname("r")
        ax.plot(r, hatrho2, linewidth = linewidth, linestyle = linestyle[i],
            label = r"${}$".format(name))

    ax.set_xlabel(r"$r$")
    ax.set_ylabel(r"$\hat{\rho}[\varphi_r]^2$")
    ax.margins(x = 0)
    ax.tick_params(right = True)
    ax.legend(frameon = False, borderaxespad = 0)

def plot_hatrho2_P1(ax, r):
    for pm, i in zip(phi_module_list, range(len(phi_module_list))):
        hatrho2_P1 = [mpmath.power(pm.phi(ri).hatrho() / P1.phi(ri).hatrho(), 2)
            for ri in r]
        name = pm.phi.latexname("r")
        ax.plot(r, hatrho2_P1, linewidth = linewidth, linestyle = linestyle[i],
            label = r"${}$".format(name))

    ax.set_xlabel(r"$r$")
    ax.set_ylabel(r"$\hat{{\rho}}[\varphi_r]^2/\hat{{\rho}}[{}]^2$".format(P1.phi.latexname("r")))
    ax.margins(x = 0)
    ax.tick_params(right = True)
    ax.legend(frameon = False, borderaxespad = 0)

matplotlib.rcParams["mathtext.fontset"] = "cm"
linewidth = 1
font_size = 10
plt.rc("font", size = font_size)
width = 5
height = 2
plt.figure(figsize = (width, height))
linestyle = ["solid", "dotted", "dashed", "dashdot"]
nrows = 1
ncols = 2
phi_module_list = [Phi, exp, P2, P1]

rmin = mpmath.power(10, -4)
rmax = mpmath.mpf("0.5")
r = [mpmath.power(10, x) for x in
    mpmath.linspace(mpmath.log10(rmin), mpmath.log10(rmax), N)]
ax = plt.subplot2grid((nrows, ncols), (0, 0), colspan = 1)
ax.set_xscale("log")
ax.set_yscale("log")
plot_hatrho2(ax, r)

rmin = mpmath.power(10, -10)
rmax = mpmath.mpf("0.99")
r = [mpmath.power(10, x) for x in
    mpmath.linspace(mpmath.log10(rmin), mpmath.log10(rmax), N)]
ax = plt.subplot2grid((nrows, ncols), (0, 1), colspan = 1)
ax.set_xlim(float(rmin), 1)
ax.set_xscale("log")
ax.set_xticks([1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1])
ax.get_xaxis().get_major_formatter().labelOnlyBase = False
plot_hatrho2_P1(ax, r)

plt.tight_layout(pad = 0, w_pad = 1)
f = __file__.replace(".py", ".pdf")
plt.savefig(f)
print(f)

