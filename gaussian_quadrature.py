import mpmath

#Compute alpha and beta of three-term recurrence relation of
#monic orthogonal polynomials induced by inner product
#<F, G> = int_a^b F(u) G(u) dW(u).
#Discretized Stieltjes procedure approximates the inner product as
#<F, G> ~ sum_{i = 1}^MDS cDSi * F(uDSi) * G(uDSi).
#alpha_m, 0 <= m <= M - 1
#beta_0 = int_a^b dW(u)
#beta_m, 1 <= m <= M - 1
def discretized_stieltjes(M, uDS, cDS):
    MDS = len(uDS)
    PP = M * [mpmath.mpf(0)]
    uPP = M * [mpmath.mpf(0)]
    alpha = M * [mpmath.mpf(0)]
    beta = M * [mpmath.mpf(0)]
    P0 =  MDS * [mpmath.mpf(1)]
    PP[0] = mpmath.fsum(cDS)
    uPP[0] = mpmath.fdot(cDS, uDS)
    alpha[0] = uPP[0] / PP[0]
    beta[0] = PP[0]

    if M == 1:
        return alpha, beta

    P1 = [u - alpha[0] for u in uDS]
    PP[1] = mpmath.fsum([c * P * P for c, P in zip(cDS, P1)])
    uPP[1] = mpmath.fsum([c * u * P * P for u, c, P in zip(uDS, cDS, P1)])
    alpha[1] = uPP[1] / PP[1]
    beta[1] = PP[1] / PP[0]

    if M == 2:
        return alpha, beta

    for m in range(2, M):
        P2 = [(ui - alpha[m - 1]) * P1i - beta[m - 1] * P0i
            for ui, P0i, P1i in zip(uDS, P0, P1)]
        PP[m] = mpmath.fsum([c * P * P for c, P in zip(cDS, P2)])
        uPP[m] = mpmath.fsum([c * u * P * P for u, c, P in zip(uDS, cDS, P2)])
        alpha[m] = uPP[m] / PP[m]
        beta[m] = PP[m] / PP[m - 1]
        P0 = P1.copy()
        P1 = P2.copy()

    return alpha, beta

def golub_welsch(alpha, beta):
    J = mpmath.zeros(len(alpha))

    for i in range(len(alpha)):
        J[i, i] = alpha[i]

    for i in range(len(alpha) - 1):
        sbeta = mpmath.sqrt(beta[i + 1])
        J[i + 1, i] = sbeta
        J[i, i + 1] = sbeta

    u, V = mpmath.eigsy(J)
    c = [beta[0] * V[0, i] * V[0, i] for i in range(len(alpha))]

    return u, c

def gauss_legendre_alpha_beta(M):
    alpha = M * [mpmath.mpf(0)]
    beta = M * [mpmath.mpf(0)]
    beta[0] = mpmath.mpf(2)

    for i in range(M - 1):
        k = mpmath.fadd(i, 1)
        k2 = k * k
        beta[i + 1] = k2 / (4 * k2 - 1)

    return alpha, beta

def gauss_legendre_quadrature_golub_welsch(M):
    alpha, beta = gauss_legendre_alpha_beta(M)
    return golub_welsch(alpha, beta)

def mpmath_quadrature_exponent(M):
    m = M
    exponent = 0

    while m % 2 == 0:
        m /= 2
        exponent += 1

    return exponent if m == 3 else -1

#M = 3 * 2^exponent
def gauss_legendre_quadrature_mpmath(M):
    exponent = mpmath_quadrature_exponent(M)

    if exponent == -1:
        raise ValueError("For mpmath, the number of nodes of the"
            " Gauss-Legendre quadrature needs to be factored as 3 * 2^m")

    nodes = mpmath.calculus.quadrature.GaussLegendre(mpmath.mp).calc_nodes(
        exponent + 1, mpmath.mp.prec)

    assert len(nodes) == M

    sorted_nodes = sorted(nodes)

    u = [u for u, c in sorted_nodes]
    c = [c for u, c in sorted_nodes]

    return u, c

