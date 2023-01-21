import mpmath
import gaussian_quadrature as gq
import util

#f(x) = sum_{i=1}^M ci * exp(-ti * x)
class ExponentialSum:
    def __init__(self, t, c):
        self.t = t.copy()
        self.c = c.copy()
        self.mct = [- c * t   for t, c in zip(self.t, self.c)]
        self.ctt = [c * t * t for t, c in zip(self.t, self.c)]

    def value(self, x):
        return mpmath.fdot(self.c,   [mpmath.exp(- t * x) for t in self.t])

    def deriv1(self, x):
        return mpmath.fdot(self.mct, [mpmath.exp(- t * x) for t in self.t])

    def deriv2(self, x):
        return mpmath.fdot(self.ctt, [mpmath.exp(- t * x) for t in self.t])

#int_{-1}^1 F(u) dW(b*phi(u))
#= int_{-1}^1 F(u) b*phi'(u) w(b*phi(u)) du
#~ sum_i cDSi F(uDSi)
def uDS_cDS_gauss_legendre(f, phi, uGL, cGL):
    uDS = uGL.copy()
    cDS = [c * f.b * phi.deriv1(u) * f.w(f.b * phi.value(u))
        for u, c in zip(uGL, cGL)]
    return uDS, cDS

def EMh_quadratic_form(f, M, h):
    F = mpmath.matrix(M + 1)
    Fij = [f.value(i * h) for i in range(2 * M + 1)]

    for i in range(M + 1):
        for j in range(M + 1):
            F[i, j] = Fij[i + j]

    u = [mpmath.mpf(1) if n % 2 == 0 else mpmath.mpf(-1) for n in range(M + 1)]

    x = mpmath.cholesky_solve(F, u)

    return 1 / mpmath.fdot(u, x)

def EMh_determinant(f, M, h):
    F = mpmath.matrix(M + 1)
    F2 = mpmath.matrix(M)
    Fi = [f.value(i * h) for i in range(2 * M + 1)]
    F2i = [Fi[i] + 2 * Fi[i+1] + Fi[i+2]  for i in range(2 * M - 1)]

    for i in range(M + 1):
        for j in range(M + 1):
            F[i, j] = Fi[i + j]

    for i in range(M):
        for j in range(M):
            F2[i, j] = F2i[i + j]

    L = mpmath.cholesky(F)
    L2 = mpmath.cholesky(F2)
    detF = mpmath.power(mpmath.fprod([L[i, i] for i in range(M + 1)]), 2)
    detF2 = mpmath.power(mpmath.fprod([L2[i, i] for i in range(M)]), 2)

    return detF / detF2

def EMh_f0_bound(a, b, M, h):
    s = mpmath.sqrt((1 + mpmath.exp(-a*h)) / (1 + mpmath.exp(-b*h)))
    return mpmath.power(s + 1, 2) * mpmath.power((s+1)/(s-1), -2*M)

#Return x satisying f(x) = y by Newton's method with the initial value x0.
def solve_y_newton(f, y, x0, print_log = False):
    #Newton's method
    #g(x) = f(x) - y
    #g'(x) = f'(x)
    #x <- x - g(x)/g'(x)
    x = x0
    flag = False
    eps = mpmath.power(2, -mpmath.mp.prec * mpmath.mpf("0.8"))
    for i in range(10000):
        dx = - (f.value(x) - y) / f.deriv1(x)
        if print_log:
            print("newton", i + 1, "dx", dx)
        x += dx
        if flag:
            return x
        elif mpmath.fabs(dx / x) <= eps:
            flag = True

    util.error("Newton's method did not converged")

def tmean(f):
    return - f.deriv1(mpmath.mpf(0)) / f.value(mpmath.mpf(0))

#Return x satisfying f(x) = y.
#The intial value of Newton's method is determined by the bisection method.
def solve_y_bisection_newton(f, y, print_log = False):
    #f(x) >= f(0) * exp(-tmean * x)
    #f(0) * exp(-tmean * x0) = y
    f0 = f.value(mpmath.mpf(0))
    x0 = mpmath.log(f0 / y) / tmean(f)
    flag = False

    for i in range(10000):
        if print_log:
            print("doubling i", i, "x0", x0)
        x1 = mpmath.ldexp(x0, 1)
        if f.value(x1) < y:
            flag = True
            break

        x0 = x1

    if flag == False:
        util.error("x0 satisfying f(x0)>=y and f(2*x0)<y could not be found")

    #Bisection method
    for i in range(10):
        xm = mpmath.ldexp(x0 + x1, -1)
        fxm = f.value(xm)
        if fxm < y:
            x1 = xm
        else:
            x0 = xm

        if print_log:
            print("bisection", i + 1, "x0", x0, "x1", x1) 

    return solve_y_newton(f, y, x0, print_log)

#f(0) - sum_{nu=1}^M c_{h,nu} <= A * B^-M * f(0)
#F = r * exp((1-r)*hr)
#A = (sqrt(F)+1)^2
#B = ((sqrt(F)+1)/(sqrt(F)-1))^2
#return hr, F, A, B
def hr_factors(r, print_log = False):
    #(1-r)*exp(-hr) + exp(-(1-r)*hr) = r
    one = mpmath.mpf(1)
    onemr = 1 - r
    t = [one, onemr]
    c = [onemr, one]
    f = ExponentialSum(t, c)
    twomr = 2 - r
    h0 = mpmath.ldexp(twomr, -1) * mpmath.log(twomr / r) / onemr
    hr = solve_y_newton(f, r, h0, print_log)
    F = r * mpmath.exp(onemr * hr)
    sF = mpmath.sqrt(F)
    p = sF + 1
    m = sF - 1
    A = p * p
    B = A / (m * m)
    return hr, F, A, B

def hr(r, print_log = False):
    _hr, F, A, B = hr_factors(r, print_log)
    return _hr

def best_Fst(f, t, c, x):
    fM = ExponentialSum(t, c)
    return [f.deriv1(xi) - fM.deriv1(xi) for xi in x]

def best_Feq(f, t, c, x):
    fM = ExponentialSum(t, c)
    E = lambda x: f.value(x) - fM.value(x)
    Ex = [E(mpmath.mpf(0))] + [E(xi) for xi in x]
    return [Ex[i] + Ex[i + 1] for i in range(len(x))]

def best_epseq(f, t, c, x):
    fM = ExponentialSum(t, c)
    E = lambda x: f.value(x) - fM.value(x)
    Ex = [E(mpmath.mpf(0))] + [E(xi) for xi in x]
    return mpmath.norm(
        [(Ex[i] + Ex[i + 1]) / Ex[i + 1] for i in range(len(x))], mpmath.inf)

def best_epsst(f, t, c, x):
    fM = ExponentialSum(t, c)
    E = lambda x: f.value(x) - fM.value(x)
    dE = lambda x: f.deriv1(x) - fM.deriv1(x)
    return mpmath.norm([xi * dE(xi) / E(xi) for xi in x], mpmath.inf)

def best_F(f, t, c, x):
    return best_Fst(f, t, c, x) + best_Feq(f, t, c, x)

def best_dFeq_dtc(t, c, x):
    M = len(t)
    dFeq_dtc = mpmath.zeros(2 * M)
    x0 = [mpmath.mpf(0)] + x

    for nu in range(M):
        emtx = [mpmath.exp(- t[nu] * xi) for xi in x0]

        for i in range(2 * M):
            dFeq_dtc[i, nu] = c[nu] * (x0[i] * emtx[i] + x0[i + 1] * emtx[i + 1])
            dFeq_dtc[i, nu + M] = - emtx[i] - emtx[i + 1]

    return dFeq_dtc

def best_dFeq_dx(f, t, c, x):
    M = len(t)
    dFeq_dx = mpmath.zeros(2 * M)
    fM = ExponentialSum(t, c)

    for i in range(2 * M):
        dFeq_dx[i, i] = f.deriv1(x[i]) - fM.deriv1(x[i])
        if i != 0:
            dFeq_dx[i, i - 1] = f.deriv1(x[i - 1]) - fM.deriv1(x[i - 1])

    return dFeq_dx

def best_dFst_dtc(t, c, x):
    M = len(t)
    dFst_dtc = mpmath.zeros(2 * M)

    for i in range(2 * M):
        for nu in range(M):
            emtx = mpmath.exp(- t[nu] * x[i])
            dFst_dtc[i, nu] = c[nu] * (1 - t[nu] * x[i]) * emtx
            dFst_dtc[i, nu + M] = t[nu] * emtx

    return dFst_dtc

def best_dFst_dx(f, t, c, x):
    M = len(t)
    dFst_dx = mpmath.zeros(2 * M)
    fM = ExponentialSum(t, c)

    for i in range(2 * M):
        dFst_dx[i, i] = f.deriv2(x[i]) - fM.deriv2(x[i])

    return dFst_dx

def best_dF_dtcx(f, t, c, x):
    M = len(t)
    dF_dtcx = mpmath.zeros(4 * M)
    M2 = 2 * M
    M4 = 4 * M

    dF_dtcx[0 : M2, 0 : M2]   = best_dFst_dtc(t, c, x)
    dF_dtcx[0 : M2, M2 : M4]  = best_dFst_dx(f, t, c, x)
    dF_dtcx[M2 : M4, 0 : M2]  = best_dFeq_dtc(t, c, x)
    dF_dtcx[M2 : M4, M2 : M4] = best_dFeq_dx(f, t, c, x)

    return dF_dtcx

def best_summary(title, f, t, c, x, print_log):
    E0 = f.value(mpmath.mpf(0)) - mpmath.fsum(c)
    epseq = best_epseq(f, t, c, x)
    epsst = best_epsst(f, t, c, x)
    s = ""

    if print_log >= 1:
        s = title + " E0 " + mpmath.nstr(E0)
        s += " max|(Ex{i-1}+Exi)/Exi| " + mpmath.nstr(epseq)
        s += " max|xi*E'xi/Exi| " + mpmath.nstr(epsst)

    if print_log >= 2:
        s += " t " + mpmath.nstr(t)
        s += " c " + mpmath.nstr(c)
        s += " x " + mpmath.nstr(x)

    if print_log >= 1:
        print(s)

    return E0, epseq, epsst

def best_newton(f, t0, c0, x0, max_iter, print_log):
    t = t0.copy()
    c = c0.copy()
    x = x0.copy()
    M = len(t)

    for j in range(max_iter + 1):
        if print_log >= 1:
            best_summary("best_newton " + str(j), f, t, c, x, print_log)

        if j == max_iter:
            return t, c, x

        mF = [-F for F in best_F(f, t, c, x)]
        dF_dtcx = best_dF_dtcx(f, t, c, x)
        dtcx = mpmath.lu_solve(dF_dtcx, mF)

        t = [t[nu] + dtcx[nu] for nu in range(M)]
        c = [c[nu] + dtcx[nu + M] for nu in range(M)]
        x = [x[i] + dtcx[i + 2 * M] for i in range(2 * M)]

#int_{e^{-bh}}^{e^{-ah}} F(y) (1+y) dWh(y)
#= int_a^b F(e^{-ht}) (1+e^{-ht}) dW(t)
#= int_{-1}^1 F(e^{-h*b*phi(u)}) * (1+e^{-h*b*phi(u)}) b*phi'(u) w(b*phi(u)) du
#~ sum_i dDSi F(yDSi)
def init_remez_yDS_dDS_gauss_legendre(f, h, phi, uGL, cGL):
    tDS = [f.b * phi.value(u) for u in uGL]
    yDS = [mpmath.exp(- h * t) for t in tDS]
    dDS = [c * (1 + y) * f.b * phi.deriv1(u) * f.w(t)
        for u, c, t, y in zip(uGL, cGL, tDS, yDS)]
    return yDS, dDS

def init_remez_gauss(alpha, beta, h):
    M = len(alpha)
    y, d = gq.golub_welsch(alpha, beta)
    t = [mpmath.log(1 / yi) / h for yi in y]
    c = [di / (1 + yi) for yi, di in zip(y, d)]
    x = mpmath.linspace(h, 2 * M * h, 2 * M)

    t.reverse()
    c.reverse()

    return t, c, x

def best_remez_update_x(f, t, c, x, max_scale, armijo, eps, print_log):
    fM = ExponentialSum(t, c)
    E  = lambda x: f.value(x)  - fM.value(x)
    x0 = x.copy()
    epsst = mpmath.mpf(0)
    nchanged = 0
    nnewton = 0

    for i in range(len(x)):
        E0 = E(x0[i])
        f1 = f.deriv1(x0[i])
        fM1 = fM.deriv1(x0[i])
        E1 = f1 - fM1
        E2 = f.deriv2(x0[i]) - fM.deriv2(x0[i])
        epsi = mpmath.fabs(x0[i] * E1 / E0)

        if epsst < epsi:
            epsst = epsi

        if epsi <= eps:
            continue

        if i == 0:
            l1 = x0[0]
            l2 = x0[1] - x0[0]
        elif i == len(x0) - 1:
            l1 = x0[i] - x0[i - 1]
            l2 = l1
        else:
            l1 = x0[i] - x0[i - 1]
            l2 = x0[i + 1] - x0[i]

        if l1 < l2:
            dx = mpmath.sign(E0 * E1) * max_scale * l1
        else:
            dx = mpmath.sign(E0 * E1) * max_scale * l2

        if E0 * E2 < 0:
            dnewton = - E1 / E2
            if  mpmath.fabs(dnewton) < mpmath.fabs(dx):
                dx = dnewton
                nnewton += 1

        if epsi > eps:
            if E0 < 0:
                dE = E(x0[i] + dx) - E0
                while dE > armijo * E1 * dx:
                    dx /= 2
                    dE = E(x0[i] + dx) - E0
            else:
                dE = E(x0[i] + dx) - E0
                while dE < armijo * E1 * dx:
                    dx /= 2
                    dE = E(x0[i] + dx) - E0
        else:
            if E0 * E2 >= 0:
                raise ValueError("best_remez_update_x"
                    " E'(xi) ~ 0 and E(xi)*E''(xi) >= 0")


        x[i] = x0[i] + dx
        nchanged += 1

    if print_log >= 1:
        best_summary("best_remez_update_x newton/changed/all "
            + str(nnewton) +  "/" + str(nchanged) + "/" + str(len(x)),
            f, t, c, x, print_log)

    return nchanged

def best_remez_update_t_c(f, t, c, x, eps, max_iter, print_log):
    for i in range(max_iter + 1):
        E0, epseq, epsst = best_summary("best_remez_update_t_c " + str(i),
            f, t, c, x, print_log)

        if epseq <= eps:
            return

        if i == max_iter:
            raise ValueError("best_remez_update_t_c"
                " did not converge for " +  str(max_iter) + " iteration")

        #Newton method
        mFeq = [-Feq for Feq in best_Feq(f, t, c, x)]
        dFeq_dtc = best_dFeq_dtc(t, c, x)
        dtc = mpmath.lu_solve(dFeq_dtc, mFeq)
        M = len(t)

        for nu in range(M):
            t[nu] += dtc[nu]
            c[nu] += dtc[nu + M]

def best_remez(f, t0, c0, x0, max_scale, armijo, eps, max_iter, print_log):
    t = t0.copy()
    c = c0.copy()
    x = x0.copy()
    M = len(t)

    for i in range(max_iter + 1):
        E0, epseq, epsst = best_summary("best_remez " + str(i), f, t, c, x, print_log)

        if epseq <= eps and epsst <= eps:
            return t, c, x, i

        if i == max_iter:
            raise ValueError("best_remez did not converge for "
                + str(max_iter) + " iteration")

        nchanged = best_remez_update_x(f, t, c, x, max_scale, armijo, eps,
            print_log - 1)

        best_remez_update_t_c(f, t, c, x, eps, max_iter, print_log - 1)

def find_max_abs_error(E, dE, x, cutoff = mpmath.mpf("0.01")):
    Ex = [E(xi) for xi in x]
    max_abs_Ex = mpmath.norm(Ex, mpmath.inf)
    abs_E_cutoff = max_abs_Ex * cutoff
    max_abs_E = max_abs_Ex

    for i in range(len(x) - 2):
        if Ex[i + 1] > 0 and Ex[i + 1] >= abs_E_cutoff:
            if Ex[i] < Ex[i + 1] and Ex[i + 1] > Ex[i + 2]:
                xst = mpmath.findroot(dE, (x[i], x[i + 2]), solver = "anderson")
                Est = E(xst)
                if Est > max_abs_E:
                    max_abs_E = Est

        if Ex[i + 1] < 0 and -Ex[i + 1] >= abs_E_cutoff:
            if Ex[i] > Ex[i + 1] and Ex[i + 1] < Ex[i + 2]:
                xst = mpmath.findroot(dE, (x[i], x[i + 2]), solver = "anderson")
                Est = E(xst)
                if -Est > max_abs_E:
                    max_abs_E = -Est

    return max_abs_E

def phi_basis_M(phi, n, eps):
    M0 = int(mpmath.ceil(mpmath.fdiv(n + 1, 2)))

    #|E| < (16/pi) * ((rho^n + rho^-n)/2) * rho^-2M <= eps * rho^-n
    #log((16/pi)/eps) + log(cosh(n*log(rho))) + n*log(rho) <= 2M*log(rho) 
    lr = mpmath.log(phi.hatrho())
    s = mpmath.log(mpmath.fdiv(16, mpmath.pi) / mpmath.mpf(eps))
    s += mpmath.log(mpmath.cosh(n * lr))
    s += n * lr
    s /= 2 * lr
    M1 = int(mpmath.ceil(s))

    return max(M0, M1)

def phi_basis(phi, n, x, M):
    s = mpmath.mpf(0)
    for nu in range(1, M + 1):
        y = mpmath.fdiv(2 * nu - 1, 2 * M)
        s += mpmath.exp(- x * phi.value(mpmath.cospi(y))) * mpmath.cospi(n * y)

    s /= M
    return s

def phi_basis_list(phi, n, xlist, M):
    y = [mpmath.fdiv(2 * nu - 1, 2 * M) for nu in range(1, M + 1)]
    phicos = [phi.value(mpmath.cospi(yi)) for yi in y]
    cosn = [mpmath.cospi(n * yi) for yi in y]
    return [mpmath.fdot([mpmath.exp(- x * phicosnu) for phicosnu in phicos],
        cosn) / M for x in xlist]

def phi_basis_sigma(n, u, c):
    return mpmath.fdot(c, [mpmath.chebyt(n, ui) for ui in u])

