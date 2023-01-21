import util
import mpmath
import mpmath_util
import cmf
import matplotlib.pyplot as plt
import matplotlib
import Phi

mpmath.mp.prec = 100
N = 501
eps = mpmath.power(10, -10)

util.print_sys_version()
mpmath_util.print_version()
mpmath_util.print_prec()
print("N", N)
print("eps", eps)

def plot_scaled_basis(ax, x, nmin, nmax, phi):
    hatrho = phi.hatrho()

    for n in range(nmin, nmax + 1):
        i = n - nmin
        M = cmf.phi_basis_M(phi, n, eps)
        print("n", n, "M", M)
        s = mpmath.power(hatrho, n)
        y = [s * cmf.phi_basis(phi, n, xi, M) for xi in x]
        ax.plot(x, y, linewidth = linewidth, linestyle = linestyle[i],
            label = r"$n={}$".format(n))

    ax.set_xscale("log")
    ax.margins(x = 0)
    ax.axhline(linewidth = 0.5, color = "k")

matplotlib.rcParams["mathtext.fontset"] = "cm"
linewidth = 1
linestyle = ["solid", "dotted", "dashed", "dashdot"]
font_size = 10
plt.rc("font", size = font_size)
width = 5
height = 6
plt.figure(figsize = (width, height))

log2r = [-1, -1, -10, -10]
nmin = [1, 5, 1, 5]
nmax = [4, 8, 4, 8]
nrows = len(log2r)
x = [mpmath.power(10, yi) for yi in mpmath.linspace(-2, 4, N)]
ylim = [0.8, 0.5, 0.8, 0.5]
legend_loc = ["lower right", "lower right", "upper left", "upper left"]

for i in range(nrows):
    print("r 2^" + str(log2r[i]))
    r = mpmath.ldexp(1, log2r[i])
    phi = Phi.phi(r)
    ax = plt.subplot2grid((nrows, 1), (i, 0), colspan = 1)
    plot_scaled_basis(ax, x, nmin[i], nmax[i], phi)
    ax.legend(ncol = 2, loc = legend_loc[i], frameon = False,
        columnspacing = 1, borderaxespad = 0)
    ax.text(0.93, 0.85, r"$r=2^{{{}}}$".format(log2r[i]),
         horizontalalignment='center',
         verticalalignment='center', transform=ax.transAxes)
    xlab = (i == nrows - 1)
    ax.tick_params(labelbottom = xlab)
    ax.set_yticks([-0.5, 0, 0.5])
    ax.set_ylim(-ylim[i], ylim[i])
    if xlab:
        ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$q^{-n} \tilde{\Phi}_{r,n}(x)$")

plt.tight_layout(pad = 0)
f = __file__.replace(".py", ".pdf")
plt.savefig(f)
print(f)

