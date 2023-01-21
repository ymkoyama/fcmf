import util
import matplotlib.pyplot as plt
import matplotlib

def dirname(param, a, b, prec):
    return "inverse_power" + param + "_a" + a + "_b" + b \
        + "_EMh_qf_prec" + str(prec)
param = ["0.5", "1", "2"]
log2a = [-1,    -10]
a = [str(2**e) for e in log2a]
b = "1"
prec  = [248,   184]
ymin  = [1e-50, 1e-30]
ymax  = [1e1,   1e1]
Mmin  = [1,     1]
Mmax  = [17,    17]

util.print_sys_version()

def plot_EMh_f0(ax, param, a, b, prec, Mmin, Mmax):
    d = dirname(param, a, b, prec)
    hab_b = util.read_float(d + "/hab_b.dat")[0]
    for M in range(Mmin, Mmax + 1):
        suffix = "_M" + f"{M:04}" + ".dat"
        h, EMh_f0 = util.read_float2(d + "/EMh_f0" + suffix)
        hb, EMh_f0b = util.read_float2(d + "/EMh_f0_bound" + suffix)
        ax.plot(h, EMh_f0, linewidth = linewidth, linestyle = "solid",
            color = "C0")
        ax.plot(hb, EMh_f0b, linewidth = linewidth, linestyle = "dotted",
            color = "C1")
    ax.axvline(x = hab_b, linewidth = 0.5, linestyle = "solid", color = "k")

    ax.margins(x = 0)
    ax.set_xscale("log")
    ax.set_yscale("log")

matplotlib.rcParams["mathtext.fontset"] = "cm"
linewidth = 1
font_size = 10
plt.rc("font", size = font_size)
width = 5
height = 5
plt.figure(figsize = (width, height))
nrows = len(param)
ncols = len(log2a)

for i in range(len(param)):
    for j in range(len(log2a)):
        ax = plt.subplot2grid((nrows, ncols), (i, j), colspan = 1)
        ax.set_title(r"$\eta={},a=2^{{{}}},b={}$".format(param[i], log2a[j], b))
        ax.set_ylim((ymin[j], ymax[j]))
        plot_EMh_f0(ax, param[i], a[j], b, prec[j], Mmin[j], Mmax[j])
        xlab = (i == len(param) - 1)
        ax.tick_params(labelbottom = xlab)
        if xlab:
            ax.set_xlabel(r"$h$")
        if j == 0:
            ax.set_ylabel(r"$E_{M,h}/f(0)$")

plt.tight_layout(pad = 0)
f = __file__.replace(".py", ".pdf")
plt.savefig(f)
print(f)

