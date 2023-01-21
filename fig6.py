import util
import mpmath
import mpmath_util
import matplotlib.pyplot as plt
import matplotlib
import importlib.util
import mpmath_util

def filename(param, a, b, MDS, prec, phi, name):
    return "inverse_power" + param + "_a" + a + "_b" + b \
        + "_best_init_MDS" + str(MDS) + "_prec" + str(prec) \
        + "/best_init_" + phi + "_" + name + "_M" + f"{M:04}" + ".dat"
param = ["0.5", "1", "2"]
log2a = [-1,  -10]
prec  = [248, 184]
xmin  = [10,  10]
ymin  = [-85, -65]
ymax  = [10,  10]
M = 17
MDS_list = []
MDS_list.append([24, 48, 96, 192, 384, 768])
MDS_list.append([24, 48, 96, 192, 384, 768])
phi_list = ["Phi", "exp", "P2", "P1", "R0_1"]
a = [str(2**e) for e in log2a]
b = "1"

util.print_sys_version()
mpmath_util.print_version()

def plot_max_rel_error(ax, x, y_list, phi_list):
    for i in range(len(y_list)):
        phi_module = importlib.import_module(phi_list[i])
        log10y = [mpmath.log10(yij) for yij in y_list[i]]
        ax.scatter(x, log10y, color = color[i], marker = marker[i],
            label = r"${}$".format(phi_module.phi.latexname("a/b")))
    ax.set_xscale("log")

matplotlib.rcParams["mathtext.fontset"] = "cm"
linewidth = 1
linestyle = ["dashdot", "dashed", "dotted", "solid"]
color = ["C1", "C2", "C3", "C4", "C5"]
marker = ["x", "1", "2", "3", "4"]
font_size = 10
plt.rc("font", size = font_size)
width = 5
height = 3.5
plt.figure(figsize = (width, height))

nrows = len(param)
ncols = len(a)

for i in range(nrows):
    for j in range(ncols):
        mpmath.mp.prec = prec[j]
        mpmath_util.print_prec()
        ax = plt.subplot2grid((nrows, ncols), (i, j), colspan = 1)
        y = []
        for phi in phi_list:
            tc = []
            for MDS in MDS_list[j]:
                t = mpmath_util.read_vector(
                filename(param[i], a[j], b, MDS, prec[j], phi, "t"))
                c = mpmath_util.read_vector(
                    filename(param[i], a[j], b, MDS, prec[j], phi, "c"))
                tc.append(t + c)

            y.append([mpmath_util.max_rel_error(tc[i + 1], tc[i])
                for i in range(len(MDS_list[j]) - 1)])

        print("MRE_h", "param", param[i], "a", a[j], "b", b, "M", M)
        print("row:phi", " ".join(phi_list))
        print("col:MDS", " ".join([str(s) for s in MDS_list[j][:-1]]))
        for yi in y:
            print(" ".join([mpmath.nstr(yij).ljust(15) for yij in yi]))

        ax.set_title(r"$\eta={},a=2^{{{}}},b=1,M={}$".format(
            param[i], log2a[j], M))
        xlab = (i == nrows - 1)
        ax.tick_params(labelbottom = xlab)
        ax.set_xlim(xmin[j], MDS_list[j][-1])
        ax.set_ylim(ymin[j], ymax[j])

        if i == nrows - 1:
            ax.set_xlabel(r"$M_\mathrm{DS}$")
        if j == 0:
            ax.set_ylabel(r"$\log_{10}(\mathrm{MRE}_h)$")

        plot_max_rel_error(ax, MDS_list[j][:-1], y, phi_list)

handles, labels = ax.get_legend_handles_labels()
leg = plt.legend(handles, labels, loc = "upper right",
    bbox_to_anchor = (1, 5.15), ncol = 5, columnspacing = 1, labelspacing = 0)
leg.set_in_layout(False)
plt.tight_layout(pad = 0, h_pad = 0.3, rect = (0, 0, 1, 0.9))
f = __file__.replace(".py", ".pdf")
plt.savefig(f)
print(f)

