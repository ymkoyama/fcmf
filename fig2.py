import util
import mpmath
import mpmath_util
import matplotlib.pyplot as plt
import matplotlib
import Phi

mpmath.mp.prec = 100
N = 201

util.print_sys_version()
mpmath_util.print_version()
mpmath_util.print_prec()
print("N", N)

def plot_phi(ax, u, phi_class, log2r):
    for i in range(len(log2r)):
        r = mpmath.ldexp(1, log2r[i])
        phi = phi_class.phi(r)
        y = [phi.value(ui) for ui in u]
        ax.plot(u, y, linewidth = linewidth, linestyle = linestyle[i],
            label = r"$r=2^{{{}}}$".format(log2r[i]))

    name = phi_class.phi.latexname("r")
    ax.set_xlabel(r"$u$")
    ax.set_ylabel(r"${}(u)$".format(name))
    ax.tick_params(right = True)

def plot_phi_deriv1(ax, u, phi_class, log2r):
    for i in range(len(log2r)):
        r = mpmath.ldexp(1, log2r[i])
        phi = phi_class.phi(r)
        y = [phi.deriv1(ui) for ui in u]
        ax.plot(u, y, linewidth = linewidth, linestyle = linestyle[i],
            label = r"$r=2^{{{}}}$".format(log2r[i]))

    name = phi_class.phi.latexname("r")
    ax.set_xlabel(r"$u$")
    ax.set_ylabel(r"${}'(u)$".format(name))
    ax.tick_params(right = True)

matplotlib.rcParams["mathtext.fontset"] = "cm"
linewidth = 1
linestyle = ["solid", "dotted", "dashed", "dashdot"]
font_size = 10
plt.rc("font", size = font_size)
width = 5
height = 2.4
plt.figure(figsize = (width, height))
nrows = 1
ncols = 2
log2r = [-1, -2, -5, -10]
u = mpmath.linspace(-1, 1, N)

ax = plt.subplot2grid((nrows, ncols), (0, 0), colspan = 1)
ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
ax.margins(x = 0)
plot_phi(ax, u, Phi, log2r)

ax = plt.subplot2grid((nrows, ncols), (0, 1), colspan = 1)
ax.margins(x = 0)
ax.set_ylim(-0.1, 2.1)
plot_phi_deriv1(ax, u, Phi, log2r)

handles, labels = ax.get_legend_handles_labels()
leg = plt.legend(handles, labels, loc = "upper right",
    bbox_to_anchor = (1, 1.24), ncol = 4)
leg.set_in_layout(False)
plt.tight_layout(pad = 0, rect = (0, 0, 1, 0.87))
f = __file__.replace(".py", ".pdf")
plt.savefig(f)
print(f)

