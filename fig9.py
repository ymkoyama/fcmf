import util
import mpmath
import mpmath_util
import matplotlib.pyplot as plt
import matplotlib
import cmf
import Phi

def file_best(param, a, b, MDS, prec, M):
    f = "inverse_power" + param + "_a" + a + "_b" + b
    f += "_best_MDS" + str(MDS) + "_prec" + str(prec)
    f += "/best_P2_E_M" + f"{M:04}" + ".dat"
    return f
def files_gauss(param, a, b, MDS, prec, M):
    d = "inverse_power" + param + "_a" + a + "_b" + b
    d += "_gauss_MDS" + str(MDS) + "_prec" + str(prec)
    E   = d + "/gauss_Phi_E_M" + f"{M:04}" + ".dat"
    eps = d + "/gauss_Phi_basis_eps_M" + f"{M:04}" + "_nmax0136.dat"
    return E, eps
param = "1"
log2a     = ["-1", "-10"]
prec      = [248,  184]
MDS_best  = [96,   192]
MDS_gauss = [96,   1536]
xmin      = [1e-2, 1e-2]
xmax      = [1e2,  2e4]
nbasis    = [1,    5]
a = [str(2**float(e)) for e in log2a]
b = "1"
M = 8
basis_eps = mpmath.power(10, -10)
basis_step = 2
basis_nx = 71

util.print_sys_version()
mpmath_util.print_version()

def read_error(filename):
    A = mpmath_util.read_matrix(filename)
    return mpmath_util.col(A, 0), mpmath_util.col(A, 3)

def plot_error(ax, x, E, x1, E1):
    EM = E[0]
    ax.plot(x, E, linewidth = linewidth, color = "k", linestyle = linestyle[0],
        label = "best")
    ax.axhline(-EM, linewidth = linewidth, color = "k",
        linestyle = linestyle[2], label = r"$\pm E_M(0)$ for best")
    ax.axhline(EM, linewidth = linewidth, color = "k", linestyle = linestyle[2])
    ax.plot(x1, E1, linewidth = linewidth, color = "k",
        linestyle = linestyle[1],
        label = r"${}$".format(Phi.phi.latexname("a/b")))

    ax.set_ylabel(r"$E_M(x)$")
    ax.set_xscale("log")
    ax.ticklabel_format(axis = "y", style = "sci", scilimits = (0, 0),
        useMathText = True)

    ymax = 1.05 * float(mpmath.norm(E1, mpmath.inf))
    ax.set_ylim(-ymax, ymax)
    ax.axhline(linewidth = 0.5, color = "k")

def plot_basis(ax, x, a, b, M, eps, nbasis, basis_step):
    phi = Phi.phi(a / b)
    nmin = 2 * M
    y = [mpmath.mpf(0)] * len(x)
    for n in range(nmin, nmin + nbasis):
        basis_M = cmf.phi_basis_M(phi, n, basis_eps)
        for i in range(len(x)):
            y[i] += 2 * eps[n] * cmf.phi_basis(phi, n, b * x[i], basis_M)

        d, m = divmod(n - nmin, basis_step)
        if m == 0:
            label = r"$N=2M$" if d == 0 else r"$N=2M+{}$".format(n - nmin)
            ax.scatter(x, y, marker = marker[d], label = label)

matplotlib.rcParams["mathtext.fontset"] = "cm"
linewidth = 1
font_size = 10
plt.rc("font", size = font_size)
width = 5
height = 5
plt.figure(figsize = (width, height))
color = ["C0", "C1", "C2", "C3", "C4"]
marker = [".", "1", "2", "3", "4"]
linestyle = ["dashed", "solid", "dotted", "dashdot"]

nrows = len(log2a)
ncols = 1
for i in range(len(a)):
    mpmath.mp.prec = prec[i]
    mpmath_util.print_prec()
    x, E = read_error(file_best(param, a[i], b, MDS_best[i], prec[i], M))
    gE, geps = files_gauss(param, a[i], b, MDS_gauss[i], prec[i], M)
    x1, E1 = read_error(gE)
    eps1 = mpmath_util.read_vector(geps)
    ax = plt.subplot2grid((nrows, ncols), (i, 0), colspan = 1)
    plot_error(ax, x, E, x1, E1)
    ax.set_xlim(xmin[i], xmax[i])
    ax.set_title(r"$\eta = {}, a = 2^{{{}}}, b = {}, M = {}$".format(
        param, log2a[i], b, M))
    basis_x = [mpmath.power(10, x) for x in mpmath.linspace(
        mpmath.log10(xmin[i]), mpmath.log10(xmax[i]), basis_nx)]
    plot_basis(ax, basis_x, mpmath.mpf(a[i]), mpmath.mpf(b), M, eps1,
        nbasis[i], basis_step)
ax.set_xlabel(r"$x$")

handles, labels = ax.get_legend_handles_labels()
leg = plt.legend(handles, labels, loc = "upper right",
    bbox_to_anchor = (1.04, 2.8), ncol = 3,
    columnspacing = 1, labelspacing = 0)
leg.set_in_layout(False)
plt.tight_layout(pad = 0, rect = (0, 0, 1, 0.9))
f = __file__.replace(".py", ".pdf")
plt.savefig(f)
print(f)

