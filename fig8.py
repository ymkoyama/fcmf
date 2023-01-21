import util
import importlib.util
import mpmath
import mpmath_util
import matplotlib.pyplot as plt
import matplotlib

def filename(param, a, b, method, MDS, prec, phi, M):
    d = "inverse_power" + param + "_a" + a + "_b" + b + "_" + method
    d += "_MDS" + str(MDS) + "_prec" + str(prec)
    f =  method + "_" + phi + "_max_abs_E_M" + f"{M:04}" + ".dat"
    return d + "/" + f
param = ["0.5", "1", "2"]
log2a        = [-1,  -10]
prec         = [248, 184]
ymin         = [30,  0.7]
ymax         = [180, 12]
yticks_start = [50,  1]
yticks_stop  = [200, 13]
yticks_step  = [50,  2]
a = [str(2**e) for e in log2a]
b = "1"
method =        ["best",  "gauss", "gauss",  "gauss",  "gauss",   "gauss"]
phi =           ["P2",    "Phi",   "exp",    "P2",     "P1",      "R0_1"]
color =         ["C0",    "C1",    "C2",     "C3",     "C4",      "C5"]
marker =        ["+",     "x",     "1",      "2",      "3",       "4"]
linestyle =     ["solid", "solid", "dotted", "dashed", "dashdot", "dashdot"]
MDS_list = []
MDS_list.append([96,      96,      96,       96,       96,        96])
MDS_list.append([192,     1536,    1536,     1536,     1536,      1536])
Mmax = 17

util.print_sys_version()
mpmath_util.print_version()

def plot_error_ratio(ax, Er_list, method, phi, r):
    for i in range(len(Er_list)):
        M = range(1, len(Er_list[i]) + 1)
        label = "best"
        if method[i] == "gauss":
            phii = importlib.import_module(phi[i]).phi(r)
            hatrho2 = mpmath.power(phii.hatrho(), 2)
            latex = phii.latexname("a/b")
            ax.axhline(hatrho2, linewidth = linewidth, linestyle = linestyle[i],
                color = "k", label = r"$\hat{{\rho}}[{}]^2$".format(latex),
                zorder = 1)
            label = r"${}$".format(latex)

        ax.scatter(M, Er_list[i], color = color[i], marker = marker[i], label = label)

matplotlib.rcParams["mathtext.fontset"] = "cm"
linewidth = 1
font_size = 10
plt.rc("font", size = font_size)
width = 5
height = 6
plt.figure(figsize = (width, height))

nrows = len(param)
ncols = len(a)
for i in range(nrows):
    for j in range(ncols):
        mpmath.mp.prec = prec[j]
        mpmath_util.print_prec()
        Er_list = []
        for k in range(len(phi)):
            fl = [filename(param[i], a[j], b, method[k], MDS_list[j][k],
                prec[j], phi[k], M) for M in range(1, Mmax + 1)]
            E = [mpmath_util.read_vector(f)[0] for f in fl]
            Er = [E[l] / E[l + 1] for l in range(len(E) - 1)]
            Er_list.append(Er)

        print("E_M/E_{M+1}", "param", param[i], "a", a[j], "b", b)
        print("row:method_phi", " ".join([m + "_" + p for m, p in zip(method, phi)]))
        print("col:M", " ".join([str(M) for M in range(1, Mmax)]))
        for Er in Er_list:
            print(" ".join([mpmath.nstr(Eri).ljust(8) for Eri in Er]))

        ax = plt.subplot2grid((nrows, ncols), (i, j), colspan = 1)
        ax.set_ylim(ymin[j], ymax[j])
        r = mpmath.mpf(a[j]) / mpmath.mpf(b)
        plot_error_ratio(ax, Er_list, method, phi, r)

        ax.text(0.3, 0.89, r"$\eta = {}, a = 2^{{{}}}, b = {}$".format(
            param[i], log2a[j], b), transform = ax.transAxes)
        xlab = (i == nrows - 1)
        ax.tick_params(labelbottom = xlab)
        if xlab:
            ax.set_xlabel(r"$M$")
        ax.set_xticks(range(1, Mmax, 5))
        if j == 0:
            ax.set_ylabel(r"$E_M/E_{M+1}$")
        ax.set_yticks(range(yticks_start[j], yticks_stop[j], yticks_step[j]))

handles, labels = ax.get_legend_handles_labels()
leg = plt.legend(handles, labels, loc = "upper right",
    bbox_to_anchor = (1.04, 3.65), ncol = 4, columnspacing = 1, labelspacing = 0)
leg.set_in_layout(False)
plt.tight_layout(pad = 0, rect = (0, 0, 1, 0.87))
f = __file__.replace(".py", ".pdf")
plt.savefig(f)
print(f)

