import util
import mpmath

def print_version():
    print("mpmath.__version__", mpmath.__version__)

def print_prec():
    print("mpmath.mp.prec", mpmath.mp.prec)

def col(A, j):
    return [A[i, j] for i in range(A.rows)]

def read_matrix(filename):
    print("read", filename)
    x = []
    f = open(filename, "r")
    i = 1
    ncol = 0
    for line in f:
        xi = [mpmath.mpf(v) for v in line.split()]
        if i == 1:
            ncol = len(xi)
            i += 1
        else:
            if len(xi) != ncol:
                util.error("The column number is different at"
                    " line=" + str(i) + " file=" + filename)
        x.append(xi)
    f.close()
    nrow = len(x)
    A = mpmath.matrix(nrow, ncol)
    for i in range(nrow):
        for j in range(ncol):
            A[i, j] = x[i][j]
    return A

def read_vector(filename):
    A = read_matrix(filename)
    if A.cols != 1:
        util.error("The column number is not 1 file=" + filename)
    return col(A, 0)

def max_abs_error(v0, v1):
    return mpmath.norm(
        [mpmath.fabs(v1[i] - v0[i]) for i in range(len(v0))],
        mpmath.inf)

def max_rel_error(v0, v1):
    return mpmath.norm(
        [mpmath.fabs((v1[i] - v0[i]) / v0[i]) for i in range(len(v0))],
        mpmath.inf)

def fAB_col(f, A, B):
    return [f(col(A, j), col(B, j)) for j in range(A.cols)]

