import argparse
import util
import matplotlib.pyplot as plt
import matplotlib

parser = argparse.ArgumentParser(
    description =
        "Plot EM(x) with the speficied prefix file\n",
    formatter_class = argparse.RawTextHelpFormatter)
parser.add_argument("--prefix", default = "prefix.dat",
    help = "File lists prefix")
parser.add_argument("--M", type = int, default = 1,
    help = "M")
args = parser.parse_args()

if args.M <= 0:
    util.error("M >= 1")

def plot_error(ax, x, E, Emax):
    ax.plot(x, E, linewidth = linewidth, color = "k", linestyle = "solid", label = "best")
    ax.axhline(-Emax, linewidth = linewidth, color = "k", linestyle = "dotted")
    ax.axhline(Emax, linewidth = linewidth, color = "k", linestyle = "dotted")

    ax.set_ylabel(r"$E_M(x)$")
    ax.set_xscale("log")
    ax.ticklabel_format(axis = "y", style = "sci", scilimits = (0, 0),
        useMathText = True)
    s = 1.2
    ax.set_ylim(- s * Emax, s * Emax)
    ax.axhline(linewidth = 0.5, color = "k")

matplotlib.rcParams["mathtext.fontset"] = "cm"
linewidth = 1
font_size = 10
plt.rc("font", size = font_size)
width = 5
height = 7
plt.figure(figsize = (width, height))

prefix = util.read_string(args.prefix)
M = args.M
Mstr = "M" + f"{M:04}"
files = [s + "_E_" + Mstr + ".dat" for s in prefix]
emax_files = [s + "_max_abs_E_" + Mstr + ".dat" for s in prefix]
nrows = len(files)
ncols = 1

for i in range(nrows):
    x, E = util.read_float2(files[i], 0, 3)
    Emax = util.read_float(emax_files[i])[0]
    ax = plt.subplot2grid((nrows, ncols), (i, 0), colspan = 1)
    plot_error(ax, x, E, Emax)
    title = files[i] + "\nmax_abs_error " + str(Emax)
    ax.set_title(title, fontsize = 6, loc = "right")
    if i == nrows - 1:
        ax.set_xlabel(r"$x$")
    ax.margins(x = 0)

plt.tight_layout(pad = 0)
f = "E_" + Mstr + ".pdf"
plt.savefig(f)
print(f)

