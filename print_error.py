import argparse
import util
import mpmath
import mpmath_util

parser = argparse.ArgumentParser(
    description = "Print max abs error and max relative error for each column",
    formatter_class = argparse.RawTextHelpFormatter)

parser.add_argument("--file0", default = "",
    help = "Reference input file")
parser.add_argument("--file1", default = "",
    help = "Target input file")
parser.add_argument("--prec", type = int, default = 1000,
    help = "Set mpmath.mp.prec (bit)")

args = parser.parse_args()

mpmath.mp.prec = args.prec
file0 = args.file0
file1 = args.file1

util.print_sys_version()
mpmath_util.print_version()
mpmath_util.print_prec()
A0 = mpmath_util.read_matrix(file0)
A1 = mpmath_util.read_matrix(file1)
print("file0", file0, "row", A0.rows, "column", A0.cols)
print("file1", file1, "row", A1.rows, "column", A1.cols)
if A1.rows != A0.rows:
    util.error("nrow(file1) != nrow(file0)")
if A1.cols != A0.cols:
    util.error("ncol(file1) != ncol(file0)")

max_ae = mpmath_util.fAB_col(mpmath_util.max_abs_error, A0, A1)
max_re = mpmath_util.fAB_col(mpmath_util.max_rel_error, A0, A1)

n = 15
print("error_type   ", " ".join(
    [("column" + str(j+1)).ljust(n) for j in range(A1.cols)]))
print("max_abs_error", " ".join([mpmath.nstr(e).ljust(n) for e in max_ae]))
print("max_rel_error", " ".join([mpmath.nstr(e).ljust(n) for e in max_re]))

