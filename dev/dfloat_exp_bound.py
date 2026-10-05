# Derives the static error bounds DFLOAT_EXP_EPS_N of the exponential
# kernels dN__exp_core (template.inc) and chooses the tapered Horner
# schedules.  Run with: python3 dev/dfloat_exp_bound.py [N] [--search]
# (needs mpmath)
#
# All bounds are absolute bounds on magnitudes, computed in 400-bit
# arithmetic with u = 2^-53, following the operations of the kernel one
# by one.  The kernels (kernels.inc) are modelled as follows.
#
# * An expansion is described by a list of bounds on the absolute
#   values of its components.  A canonical expansion (as verified by
#   _is_canonical) with head bound X has bounds X u^i.
#
# * Order accumulation: order k of a product or sum is the TwoSum
#   accumulation of a list of terms in a fixed order; each TwoSum
#   acc = fl(acc + t) produces a residual of magnitude at most
#   u |acc|, and |acc| <= (1 + u)(|acc_old| + |t|).  The residuals of
#   order k are terms of order k + 1 (accumulated last, in the order
#   they were produced); those of the last order are dropped.
#
# * _renorm_orders (bottom-up VecSum of S_0..S_{M-1}): T_{M-1} =
#   S_{M-1}, T_i = fl(S_i + T_{i+1}); the output is z_0 = T_0 and
#   z_{i+1} = the rounding error of T_i, so |z_{i+1}| <= u |T_i| and
#   |T_i| <= (|S_i| + |T_{i+1}|)(1 + u).  Nothing is dropped.  A zero
#   head is shifted away, which can only decrease the components.
#
# * mul: exact products x_i y_j = h_ij + l_ij for i + j < M, |h_ij| <=
#   (1 + u)|x_i y_j|, |l_ij| <= u |x_i y_j|; order k collects the h of
#   order k, the l of order k - 1 and the residuals of order k - 1.
#   Dropped: the products of order >= M, the l of order M - 1 and the
#   residuals of order M - 1.  sqr: the same with symmetric terms
#   doubled.  mul_d (times a double): products h_i + l_i.
#
# * add/sub: componentwise s_i + e_i = x_i + y_i, |s_i| <= (1 + u)
#   (|x_i| + |y_i|), |e_i| <= u |s_i| (s_i = x_i + y_i exactly and
#   e_i = 0 when one of them is zero); order k collects s_k, e_{k-1}
#   and the residuals of order k - 1.  Dropped: e_{M-1} and the
#   residuals of order M - 1.  add3 (three terms): two TwoSums per
#   component, both residuals go to the next order.
#
# Subnormal intermediate quantities (below 2^-960 in magnitude) can
# make an exact transformation inexact by at most 2^-1074; a blanket
# 2^-1000 absolute term covers the at most 2^8 such operations.

import sys
from math import factorial, exp, log, ceil, log2, sin
from mpmath import mp, mpf

# The arithmetic below is done with 400-bit floating-point numbers
# (mpmath); every operation is a sum or product of positive bounds,
# so the relative rounding error is at most 2^-390 per operation,
# which the final inflation by (1 + 2^-300) absorbs.
mp.prec = 400

def F(a, b=1):
    return mpf(a) / b

u = F(1, 2**53)
ONE_PLUS_U = 1 + u

def hexfloat(s):
    return F(float.fromhex(s))

# constants of the tables (exp_tables.c)
L = [hexfloat(t) for t in ["0x1.62e42fefa4000p-1", "-0x1.8432a1b0e2634p-43",
     "0x1.f97b57a079a19p-103", "0x1.9ca62d8b62834p-158",
     "0x1.75b8baafa2be8p-212", "-0x1.1e277e48514dap-266"]]
LN2_RESID = [None, F(1, 2**102), F(1, 2**157), F(1, 2**211), F(1, 2**265), F(1, 2**321)]
TAB_ERR = F(1, 2**268)
KMAX = 1078
X0MAX = F(747)

def fr(x):
    return "2^%.2f" % log2(float(x)) if x > 0 else "0"

# ---- order accumulation ----

def accumulate(orders, M):
    """orders: list of M lists of term bounds, in the order the kernel
    accumulates them (before the residuals of the previous order,
    which come last).  Returns the order sums S (bounds) and the
    dropped residual bound of the last order.  The residuals of order
    k feed order k + 1."""
    S = []
    resid_prev = []
    for k in range(M):
        terms = list(orders[k]) + resid_prev
        acc = terms[0]
        resid = []
        for t in terms[1:]:
            acc = ONE_PLUS_U * (acc + t)
            resid.append(u * acc)
        S.append(acc)
        resid_prev = resid
    return S, sum(resid_prev, F(0))

def renorm(S):
    M = len(S)
    T = [F(0)] * M
    T[M - 1] = S[M - 1]
    for i in range(M - 2, -1, -1):
        T[i] = (S[i] + T[i + 1]) * ONE_PLUS_U
    z = [T[0]] + [u * T[i] for i in range(M - 1)]
    # monotone, so that a zero-head shift keeps the bounds valid
    for i in range(M - 1):
        assert z[i] >= z[i + 1]
    return z

def mul(x, y, M):
    """x, y: component bounds (length M).  Returns (z bounds, dropped)."""
    orders = [[] for _ in range(M)]
    for i in range(M):
        for j in range(M - i):
            orders[i + j].append(ONE_PLUS_U * x[i] * y[j])
            if i + j + 1 < M:
                orders[i + j + 1].append(u * x[i] * y[j])
    S, resid = accumulate(orders, M)
    dropped = resid
    for i in range(M):
        for j in range(M):
            if i + j >= M:
                dropped += x[i] * y[j]
            elif i + j == M - 1:
                dropped += u * x[i] * y[j]
    return renorm(S), dropped

def sqr(x, M):
    orders = [[] for _ in range(M)]
    for i in range(M):
        for j in range(i, M - i):
            f = 1 if i == j else 2
            orders[i + j].append(f * ONE_PLUS_U * x[i] * x[j])
            if i + j + 1 < M:
                orders[i + j + 1].append(f * u * x[i] * x[j])
    S, resid = accumulate(orders, M)
    dropped = resid
    for i in range(M):
        for j in range(i, M):
            f = 1 if i == j else 2
            if i + j >= M:
                dropped += f * x[i] * x[j]
            elif i + j == M - 1:
                dropped += f * u * x[i] * x[j]
    return renorm(S), dropped

def mul_d(x, y, M):
    """x: component bounds, y: a bound on the double factor."""
    orders = [[] for _ in range(M)]
    for i in range(M):
        orders[i].append(ONE_PLUS_U * x[i] * y)
        if i + 1 < M:
            orders[i + 1].append(u * x[i] * y)
    S, resid = accumulate(orders, M)
    return renorm(S), resid + u * x[M - 1] * y

def add(x, y, M, xzero=(), yzero=()):
    """xzero / yzero: indices of components known to be exactly zero."""
    orders = [[] for _ in range(M)]
    for i in range(M):
        if i in xzero or i in yzero:
            s_i = x[i] + y[i]
            e_i = F(0)
        else:
            s_i = ONE_PLUS_U * (x[i] + y[i])
            e_i = u * s_i
        orders[i].append(s_i)
        if i + 1 < M:
            orders[i + 1].append(e_i)
        else:
            dropped_e = e_i
    S, resid = accumulate(orders, M)
    return renorm(S), resid + dropped_e

def add3(x, y, z, M):
    """x + y + z: componentwise s_i + e1_i + e2_i = x_i + y_i + z_i by
    two TwoSums; order k collects s_k, e1_{k-1}, e2_{k-1} and the
    residuals of order k - 1.  Dropped: e1_{M-1}, e2_{M-1} and the
    residuals of order M - 1."""
    orders = [[] for _ in range(M)]
    for i in range(M):
        s1 = ONE_PLUS_U * (x[i] + y[i])
        e1 = u * s1
        s2 = ONE_PLUS_U * (s1 + z[i])
        e2 = u * s2
        orders[i].append(s2)
        if i + 1 < M:
            orders[i + 1].append(e1)
            orders[i + 1].append(e2)
        else:
            dropped_e = e1 + e2
    S, resid = accumulate(orders, M)
    return renorm(S), resid + dropped_e

def addn(xs, M):
    """sum of the expansions xs (the addn kernel): componentwise, the
    nx - 1 TwoSums' residuals all go to the next order."""
    orders = [[] for _ in range(M)]
    dropped_e = F(0)
    for i in range(M):
        acc = xs[0][i]
        es = []
        for x in xs[1:]:
            acc = ONE_PLUS_U * (acc + x[i])
            es.append(u * acc)
        orders[i].append(acc)
        if i + 1 < M:
            orders[i + 1].extend(es)
        else:
            dropped_e = sum(es, F(0))
    S, resid = accumulate(orders, M)
    return renorm(S), resid + dropped_e

def canonical(head, M):
    return [head * u ** i for i in range(M)]

def truncate(x, M):
    """the first M components and the magnitude of the rest"""
    return x[:M], sum(x[M:], F(0))

def const_bounds(k, M):
    # nearest expansion of 1/k!: |c_i| <= u^i (1/k!)(1 + u), 5 terms
    c = F(1, factorial(k)) * ONE_PLUS_U
    comps = [c * u ** i for i in range(5)]
    return comps[:M], sum(comps[M:], F(0)) + TAB_ERR

def poly_deriv(degs, Z):
    """|d/dz sum_m c_{degs[m]} z^m| <= this for |z| <= Z, c_k = 1/k!:
    the sensitivity of a chain to the error of its variable (the
    chains are evaluated at the computed w or v, whose errors e_w, e_v
    against s^2, w^2 enter through this)"""
    return sum((m * F(1, factorial(degs[m])) * Z ** (m - 1) for m in range(1, len(degs))), F(0)) * ONE_PLUS_U

# ---- the kernel ----

def bound(N, lev, split=2, verbose=False, m1=None):
    """m1 = (regime, param): the expm1 variant (k = 0, |x| <= 0.34,
    the constant term C = E1 + T1 E2 instead of T); regime 'a': i =
    j = 0, param = X <= 2^-17 (the bound is for |x| in [X / ratio,
    X]); 'b': i = 0, j != 0, param = X in [2^-17, 2^-9]; 'c': param =
    i >= 1 (|x| within 1/512 of i/256).  Returns (E_total, R_min)
    then.

    split = 2 (the vector kernel): lev[j], j = 0..JE, the level of
    step j of the chains E and O in w (the O chain, shorter by
    JE - JO, is padded with leading zero coefficients, which are exact
    and contribute nothing, so it is modelled with its own JO + 1
    steps at the levels lev[0..JO]); split = 4 (the scalar kernel):
    lev[j], j = 0..J4, the level of step j of the four chains E0, E1,
    O0, O1 in v = w^2, and y = T + (T s) E0 + (T s w) E1 + (T w) O0 +
    (T w^2) O1 as one five-term sum."""
    deg = {1: 3, 2: 6, 3: 8, 4: 11}[N]
    JE = (deg - 1) // 2
    JO = (deg - 2) // 2
    if split == 2:
        assert len(lev) == JE + 1 and lev[0] == N
    else:
        J4 = (JE + 1 + 1) // 2 - 1
        assert len(lev) == J4 + 1 and lev[0] == N

    if m1 is None:
        # 1. canonical input
        xc = canonical(X0MAX, N)

        # 2. reduction: |t| <= log(2)/2 + 2^-43, |h| <= |t| + k|log 2 - L0| + tail
        S0 = F(1, 2**17)
        # 3. Kt = k (L1 + ... + LN)
        Lb = [abs(L[i]) for i in range(1, N + 1)]
        Kt, e_Kt = mul_d(Lb, F(KMAX), N)
        e_ln2 = KMAX * LN2_RESID[N]
        # 4. Y = (x1, ..., x_{N-1}, 0) - Kt
        A = xc[1:] + [F(0)]
        Y, e_Y = add(A, Kt, N, xzero=(N - 1,))
        # 5. s = (s0, 0, ...) + Y (the component bounds of the weak output
        # are used as they are)
        B = [S0] + [F(0)] * (N - 1)
        s, e_add = add(B, Y, N, xzero=tuple(range(1, N)))
        E_s = e_ln2 + e_Kt + e_Y + e_add
        Smax = S0 + sum(Y, F(0))          # |s|

        # 6. T = T1 T2 (tables truncated to N terms), T s, T s^2
        T1max = F(exp(89 / 256)) * F(1001, 1000)
        T2max = F(exp(128 / 65536)) * F(1001, 1000)
        T1 = canonical(T1max * ONE_PLUS_U, 5)
        T2 = canonical(T2max * ONE_PLUS_U, 5)
        T1N, e_T1 = truncate(T1, N)
        T2N, e_T2 = truncate(T2, N)
        e_T1 += TAB_ERR
        e_T2 += TAB_ERR
        T, e_mT = mul(T1N, T2N, N)
        Tval = sum(T, F(0))
        E_T = e_mT + e_T1 * T2max + e_T2 * T1max
    else:
        # expm1: k = 0 (Kt = 0 and Y = the tail of x exactly), s = s0 + Y
        regime, param = m1
        if regime == 'a':
            X = F(param); S0 = X
            xc = canonical(X, N)
            T1max = F(1); T2max = F(1)
            T = [F(1)] + [F(0)] * (N - 1); Tval = F(1); E_T = F(0)
            C = [F(0)] * N; E_C = F(0)
            Rmin = X / (mpf(2) ** F(1, 4)) * (1 - X / 2)
        elif regime == 'b':
            X = F(param); S0 = F(1, 2**17)
            xc = canonical(X, N)
            T1max = F(1)
            E2max = X + F(1, 2**17)                       # |exp(j/65536) - 1|
            T2max = 1 + E2max
            T2 = canonical(T2max * ONE_PLUS_U, 5)
            T2N, e_T2 = truncate(T2, N); e_T2 += TAB_ERR
            T = T2N; Tval = T2max; E_T = e_T2
            # C = E1 + T1 E2 = E2 exactly (T1 = 1, E1 = 0)
            C, E_C = coef_bounds(E2max, N)
            Rmin = X / (mpf(2) ** F(1, 4)) * (1 - X / 2)
        else:
            i = param; S0 = F(1, 2**17)
            X = F(i, 256) + F(1, 512)
            xc = canonical(X, N)
            T1max = F(exp(i / 256)) * F(1001, 1000)
            T2max = F(exp(128 / 65536)) * F(1001, 1000)
            T1 = canonical(T1max * ONE_PLUS_U, 5); T2 = canonical(T2max * ONE_PLUS_U, 5)
            T1N, e_T1 = truncate(T1, N); T2N, e_T2 = truncate(T2, N)
            e_T1 += TAB_ERR; e_T2 += TAB_ERR
            T, e_mT = mul(T1N, T2N, N)
            Tval = sum(T, F(0))
            E_T = e_mT + e_T1 * T2max + e_T2 * T1max
            E1max = F(exp(i / 256) - 1) * F(1001, 1000)
            E2max = F(exp(128 / 65536) - 1) * F(1001, 1000)
            E1b, e_E1 = coef_bounds(E1max, N)
            E2b, e_E2 = coef_bounds(E2max, N)
            P12, e_m = mul(T1N, E2b, N)
            C, e_a = add(E1b, P12, N)
            E_C = e_a + e_E1 + e_m + e_T1 * E2max + T1max * e_E2
            # |expm1| is smallest at the negative end: 1 - exp(-(i/256 - 1/512))
            Rmin = 1 - mp.exp(-(F(i, 256) - F(1, 512)))
        A = xc[1:] + [F(0)]
        Y = A; e_Y = F(0)
        B = [S0] + [F(0)] * (N - 1)
        s, e_add = add(B, Y, N, xzero=tuple(range(1, N)))
        E_s = e_add
        Smax = S0 + sum(Y, F(0))
    # 7. w = s^2, T s, T w
    w, e_w = sqr(s, N)
    Wval = Smax ** 2                   # |w| <= this + e_w
    Ts, e_m = mul(T, s, N)
    E_Ts = e_m + E_T * Smax
    Tsval = sum(Ts, F(0))
    Ts2, e_m = mul(T, w, N)
    E_Ts2 = e_m + E_T * (Wval + e_w) + Tval * e_w
    Ts2val = sum(Ts2, F(0))

    # 8. the Horner chains: q = sum_j w^j c_j evaluated from the top;
    # track (component bounds, value bound, error bound) where the
    # error is against the exact polynomial in the computed w
    def chain(levels, coeff):
        J = len(levels) - 1
        M = levels[J]
        q, e_c = const_bounds(coeff(J), M)
        Eq = e_c
        for j in range(J - 1, -1, -1):
            M = levels[j]
            Mprev = levels[j + 1]
            assert M >= Mprev
            qM = q + [F(0)] * (M - Mprev)
            wM, e_wtrunc = truncate(w, M)
            prod, e_mul = mul(wM, qM, M)
            # exact w q vs computed: mul error, w truncation, q error
            Eprod = e_mul + e_wtrunc * sum(qM, F(0)) + (Wval + e_w) * Eq
            c, e_c = const_bounds(coeff(j), M)
            q, e_add = add(c, prod, M)
            Eq = e_add + e_c + Eprod
        return q, Eq
    if split == 2:
        qE, E_E = chain(lev, lambda j: 2 * j + 1)
        qO, E_O = chain(lev[:JO + 1], lambda j: 2 * j + 2)
        qE = qE + [F(0)] * (N - len(qE))
        qO = qO + [F(0)] * (N - len(qO))
        Eval = sum(qE, F(0))
        Oval = sum(qO, F(0))

        # 9. v0 = (T s) E, v1 = (T w) O, y = T + v0 + v1 (add3)
        v0, e_m = mul(Ts, qE, N)
        E_v0 = e_m + E_Ts * Eval + Tsval * E_E
        v1, e_m = mul(Ts2, qO, N)
        E_v1 = e_m + E_Ts2 * Oval + Ts2val * E_O
        y, e_a = add3(C if m1 else T, v0, v1, N)
        # y vs T_comp (1 + s E* + s^2 O*) where E*, O* are the exact
        # polynomials at the computed s (the errors E_T of T itself are
        # counted once, in E_v0, E_v1 and below); the chains are
        # evaluated at the computed w = s^2 + e_w
        E_y = e_a + E_v0 + E_v1
        E_y += e_w * (Tsval * poly_deriv(list(range(1, deg + 1, 2)), Wval + e_w)
                      + Ts2val * poly_deriv(list(range(2, deg + 1, 2)), Wval + e_w))
    else:
        # the chains in v = w^2: E0 = c1 + c5 v + c9 v^2 + ..., etc.
        v_, e_v = sqr(w, N)
        Vval = (Wval + e_w) ** 2
        Tsw, e_m = mul(Ts, w, N)
        E_Tsw = e_m + E_Ts * (Wval + e_w) + Tsval * e_w
        Tswval = sum(Tsw, F(0))
        Tw2, e_m = mul(Ts2, w, N)
        E_Tw2 = e_m + E_Ts2 * (Wval + e_w) + Ts2val * e_w
        Tw2val = sum(Tw2, F(0))
        def chain4(levels, coeff):
            # like chain, in v, with coefficient degrees 4j + offset,
            # zero (exact) coefficients beyond the degree
            J = len(levels) - 1
            # the highest nonzero coefficient (an all-zero lane is exact)
            nz = [j for j in range(J + 1) if coeff(j) <= deg]
            if not nz:
                return [F(0)] * N, F(0)
            Jn = max(nz)
            M = levels[Jn]
            q, e_c = const_bounds(coeff(Jn), M)
            Eq = e_c
            for j in range(Jn - 1, -1, -1):
                M = levels[j]
                Mprev = levels[j + 1]
                qM = q + [F(0)] * (M - Mprev)
                vM, e_vtrunc = truncate(v_, M)
                prod, e_mul = mul(vM, qM, M)
                Eprod = e_mul + e_vtrunc * sum(qM, F(0)) + (Vval + e_v) * Eq
                c, e_c = const_bounds(coeff(j), M)
                q, e_add = add(c, prod, M)
                Eq = e_add + e_c + Eprod
            return q + [F(0)] * (N - len(q)), Eq
        qE0, E_E0 = chain4(lev, lambda j: 4 * j + 1)
        qE1, E_E1 = chain4(lev, lambda j: 4 * j + 3)
        qO0, E_O0 = chain4(lev, lambda j: 4 * j + 2)
        qO1, E_O1 = chain4(lev, lambda j: 4 * j + 4)
        v0, e_m = mul(Ts, qE0, N)
        E_v0 = e_m + E_Ts * sum(qE0, F(0)) + Tsval * E_E0
        v1, e_m = mul(Tsw, qE1, N)
        E_v1 = e_m + E_Tsw * sum(qE1, F(0)) + Tswval * E_E1
        v2, e_m = mul(Ts2, qO0, N)
        E_v2 = e_m + E_Ts2 * sum(qO0, F(0)) + Ts2val * E_O0
        v3, e_m = mul(Tw2, qO1, N)
        E_v3 = e_m + E_Tw2 * sum(qO1, F(0)) + Tw2val * E_O1
        y, e_a = addn([C if m1 else T, v0, v1, v2, v3], N)
        E_y = e_a + E_v0 + E_v1 + E_v2 + E_v3
        # the chains at the computed v = w^2 + e_v, and E = E0 + w E1,
        # O = O0 + w O1 at the computed w = s^2 + e_w
        E_y += e_v * (Tsval * poly_deriv(list(range(1, deg + 1, 4)), Vval + e_v)
                      + Tswval * poly_deriv(list(range(3, deg + 1, 4)), Vval + e_v)
                      + Ts2val * poly_deriv(list(range(2, deg + 1, 4)), Vval + e_v)
                      + Tw2val * poly_deriv(list(range(4, deg + 1, 4)), Vval + e_v))
        E_y += e_w * (Tsval * poly_deriv(list(range(1, deg + 1, 2)), Wval + e_w)
                      + Ts2val * poly_deriv(list(range(2, deg + 1, 2)), Wval + e_w))
        E_E = E_E0 + E_E1
        E_O = E_O0 + E_O1

    # Taylor remainder of expm1 beyond degree deg, at |s| <= Smax
    Sf = float(Smax)
    R_T = F(Sf ** (deg + 1) / factorial(deg + 1) * exp(Sf)) * F(1001, 1000)

    # y_true = T_true exp(s_true); exp(s_true) = (1 + p* + R_T) exp(s_true - s_comp)
    eS = F(exp(Sf)) * F(1001, 1000)
    e_red = eS * F(exp(float(E_s)) - 1) * F(1001, 1000) + F(1, 2**1000)
    if m1:
        # expm1(x) = C_true + T_true (exp(s_true) - 1): |y - expm1(x)| <=
        # E_y (against C + T (p*)) + E_C + T (R_T + e_red) + E_T |exp(s) - 1|
        E_total = E_y + E_C + Tval * (R_T + e_red) + E_T * (eS - 1 + F(1, 2**1000)) + F(1, 2**1000)
        if verbose:
            print("  expm1 N = %d split %d regime %s param %s: error %s, |expm1| >= %s, rel %s" % (N, split, m1[0], m1[1], fr(E_total), fr(Rmin), fr(E_total / Rmin)))
        return E_total * (1 + F(1, 2**300)), Rmin
    # |y - T_true exp(s_true)| <= E_y (against T_comp (1 + p*)) + |T_comp| (R_T + e_red) + E_T exp(s_true)
    E_total = E_y + Tval * (R_T + e_red) + E_T * eS + F(1, 2**1000)
    ymin = F(exp(-89 / 256 - 128 / 65536)) * F(999, 1000) / eS
    rel = E_total / ymin * (1 + F(1, 2**300))

    if verbose:
        print("N = %d, deg = %d, split %d, levels %s" % (N, deg, split, lev))
        print("  |s| <= %s, reduction error %s (ln2 %s, Kt %s, Y %s, add %s)" % (fr(Smax), fr(E_s), fr(e_ln2), fr(e_Kt), fr(e_Y), fr(e_add)))
        print("  w error %s, E chain error %s, O chain error %s" % (fr(e_w), fr(E_E), fr(E_O)))
        print("  T error %s, Ts %s, Ts2 %s, v0 %s, v1 %s, y error %s, Taylor %s" % (fr(E_T), fr(E_Ts), fr(E_Ts2), fr(E_v0), fr(E_v1), fr(E_y), fr(R_T)))
        print("  total absolute %s, relative %s" % (fr(E_total), fr(rel)))
    return rel

def bound1(verbose=False):
    """The double kernel (N = 1): h = x0 - k L0 and d = h - i/256
    exactly, s = fl(fl(d - p1) - e1) with p1 + e1 = k L1 exactly,
    P = fma chain of the Taylor coefficients 1/720 ... 1 (rounded
    doubles), p = fl(s P), y = fma(T, p, T) with T the nearest double
    of exp(i/256); every operation rounds once with relative error u."""
    deg = 6
    Sd = F(1, 2**9) + F(1078) * abs(L[1])      # |d - p1| <= this
    S = Sd * ONE_PLUS_U + u * Sd                # |s_comp|
    # reduction error: two roundings, dropped k L2 + ...
    E_s = u * Sd + u * S + KMAX * LN2_RESID[1]
    # the polynomial: P_j = fma(s, P_{j+1}, c_j), c_j = rounded 1/j!
    P = F(1, factorial(deg)) * ONE_PLUS_U
    EP = u * F(1, factorial(deg))
    for j in range(deg - 1, 0, -1):
        c = F(1, factorial(j)) * ONE_PLUS_U
        val = (S * P + c) * ONE_PLUS_U
        EP = u * val + S * EP + u * F(1, factorial(j))
        P = val
    p = S * P * ONE_PLUS_U
    E_p = u * S * P + S * EP
    Sf = float(S)
    R_T = F(Sf ** (deg + 1) / factorial(deg + 1) * exp(Sf)) * F(1001, 1000)
    T1max = F(exp(89 / 256)) * F(1001, 1000)
    E_T = u * T1max
    y = T1max * (1 + p) * ONE_PLUS_U
    E_y = u * y + T1max * E_p
    eS = F(exp(Sf)) * F(1001, 1000)
    e_red = eS * F(exp(float(E_s)) - 1) * F(1001, 1000) + F(1, 2**1000)
    E_total = E_y + T1max * (R_T + e_red) + E_T * eS + F(1, 2**1000)
    ymin = F(exp(-89 / 256)) * F(999, 1000) / eS
    rel = E_total / ymin * (1 + F(1, 2**300))
    if verbose:
        print("N = 1 (double kernel): |s| <= %s, reduction error %s, p error %s, Taylor %s, y error %s" % (fr(S), fr(E_s), fr(E_p), fr(R_T), fr(E_y)))
        print("  total absolute %s, relative %s" % (fr(E_total), fr(rel)))
    return rel

def expm1_eps(N, split, lev, verbose=False):
    """EPS of expm1: the kernel's regimes for |x| <= 0.34, and
    exp(x) - 1 beyond (|exp(x) - 1| >= 0.29: the exponential's
    relative bound times exp(x) / |exp(x) - 1| <= 3.4, plus the
    subtraction's u^N |y| against |y - 1|)"""
    eps = F(0)
    ratio = mpf(2) ** F(1, 4)
    X = mpf(2) ** -17
    while X > mpf(2) ** -300:
        E, R = bound(N, lev, split, m1=('a', X))
        eps = max(eps, E / R)
        X /= ratio
    X = mpf(2) ** -17 * ratio
    while X < mpf(2) ** -9 * ratio:
        E, R = bound(N, lev, split, m1=('b', min(X, mpf(2) ** -9)))
        eps = max(eps, E / R)
        X *= ratio
    for i in range(1, 90):
        E, R = bound(N, lev, split, m1=('c', i))
        eps = max(eps, E / R)
    # beyond 0.34
    ee = bound(N, LEVELS[N] if split == 2 else LEVELS4[N], split)
    e1 = ee * F(exp(0.34)) / F(exp(0.34) - 1)
    e2 = ee * F(exp(-0.34)) / F(1 - exp(-0.34))
    # the subtraction y - 1 (an add kernel): u^N |y| roughly, against |y - 1|
    yb = canonical(F(exp(0.34)) * F(1001, 1000), N)
    one = [F(1)] + [F(0)] * (N - 1)
    d, e_sub = add(yb, one, N)
    eps = max(eps, max(e1, e2) + e_sub / F(exp(0.34) - 1))
    if verbose:
        print("N = %d split %d levels %s: EXPM1 EPS %s" % (N, split, lev, fr(eps)))
    return eps

def bound1_expm1(regime, param, verbose=False):
    """The double kernel: as bound1 with T (1 + p) replaced by fma(T, p,
    E1) (E1 = exp(i/256) - 1 rounded, k = 0 and no k L1 term)."""
    deg = 6
    if regime == 'a':
        X = F(param); S = X; E_s = F(0)
        T = F(1); E_T = F(0); E1 = F(0); E_E1 = F(0)
        Rmin = X / (mpf(2) ** F(1, 4)) * (1 - X / 2)
    else:
        i = param
        S = F(1, 512) * ONE_PLUS_U; E_s = F(0)      # s = h - i/256 exact
        T = F(exp(i / 256)) * F(1001, 1000); E_T = u * T
        E1 = F(exp(i / 256) - 1) * F(1001, 1000); E_E1 = u * E1
        Rmin = 1 - mp.exp(-(F(i, 256) - F(1, 512)))
    P = F(1, factorial(deg)) * ONE_PLUS_U
    EP = u * F(1, factorial(deg))
    for j in range(deg - 1, 0, -1):
        c = F(1, factorial(j)) * ONE_PLUS_U
        val = (S * P + c) * ONE_PLUS_U
        EP = u * val + S * EP + u * F(1, factorial(j))
        P = val
    p = S * P * ONE_PLUS_U
    E_p = u * S * P + S * EP
    Sf = float(S)
    R_T = F(Sf ** (deg + 1) / factorial(deg + 1) * exp(Sf)) * F(1001, 1000)
    y = (T * p + E1) * ONE_PLUS_U
    E_y = u * y + T * E_p + E_T * p + E_E1
    E_total = E_y + T * R_T + F(1, 2**1000)
    if verbose:
        print("  expm1 N = 1 regime %s param %s: error %s, rel %s" % (regime, param, fr(E_total), fr(E_total / Rmin)))
    return E_total * (1 + F(1, 2**300)), Rmin

def expm1_eps1(verbose=False):
    eps = F(0)
    ratio = mpf(2) ** F(1, 4)
    X = mpf(2) ** -9
    while X > mpf(2) ** -300:
        E, R = bound1_expm1('a', X)
        eps = max(eps, E / R)
        X /= ratio
    for i in range(1, 90):
        E, R = bound1_expm1('c', i)
        eps = max(eps, E / R)
    ee = bound1()
    eps = max(eps, ee * F(exp(0.34)) / F(exp(0.34) - 1) + u * F(exp(0.34)) / F(exp(0.34) - 1))
    if verbose:
        print("N = 1 (double kernel): EXPM1 EPS %s" % fr(eps))
    return eps

# ---- sine and cosine ----

PI2 = [hexfloat(t) for t in ["0x1.921fb54400000p+0", "0x1.0b4611a626331p-34",
       "0x1.1701b839a2520p-88", "0x1.27044533e63a0p-142", "0x1.05df531d89cd9p-198",
       "0x1.28a5043cc71a0p-254"]]
PI2_RESID = [None, F(1, 2**87), F(1, 2**141), F(1, 2**197), F(1, 2**253), F(1, 2**308)]
TRIG_KMAX = 2**20
TRIG_DEG = {1: (7, 8), 2: (7, 8), 3: (11, 10), 4: (13, 12)}

def bound_trig(N, split, X, regime, lev, verbose=False):
    """The absolute error of sin (and cos) of a canonical x with
    |x_0| <= X.  regime 'tiny': X < 2^-15, no reduction and no table
    (s = x exactly, SA = 0, CA = 1); 'small': 2^-15 <= X < 3/4 (k = 0:
    no reduction; the table point A = i/64 + j/2^13 has SA = sin A <=
    X + 2^-13 <= 5 X); 'general': the reduction with |k| <= 2^20 and
    SA, CA <= 1.  Returns (E_sin, E_cos)."""
    dsin, dcos = TRIG_DEG[N]
    Ks = (dsin - 1) // 2          # P: coefficients of w^0..w^Ks
    Kc = dcos // 2                # Q: coefficients of w^0..w^(Kc-1)
    X = F(X)
    xc = canonical(X, N)
    tab_trunc = lambda M: sum((u ** i for i in range(M, 5)), F(0)) * ONE_PLUS_U + TAB_ERR

    if regime == 'general':
        # h = x0 - k P0 exactly, |h| <= pi/4 + 2^-13 + 2^-30
        H = F(0.7854) + F(1, 2**13) + F(1, 2**30)
        Pb = [abs(PI2[i]) for i in range(1, N + 1)]
        Kt, e_Kt = mul_d(Pb, F(TRIG_KMAX), N)
        e_pi = TRIG_KMAX * PI2_RESID[N]
        A = [H] + xc[1:]
        T, e_T = add(A, Kt, N)
        E_red = e_pi + e_Kt + e_T           # |t_comp - t_true|
        # t is canonicalised (exactly): its tail is bounded by u^i |t|
        T = canonical(sum(T, F(0)), N)
        s = [F(1, 2**14) + F(1, 2**52)] + T[1:]
        SAmax, CAmax = F(1), F(1)
        # SA = sa cb + ca sb, CA = ca cb - sa sb from nearest 5-term
        # expansions truncated to N terms (|sa|, |sb| <= 1)
        sa = canonical(F(1) * ONE_PLUS_U, N); sb = sa; ca = sa; cb = sa
        e_tab = tab_trunc(N)
        p1, e1 = mul(sa, cb, N); p2, e2 = mul(ca, sb, N)
        SA, e3 = add(p1, p2, N)
        E_SA = e1 + e2 + e3 + 2 * e_tab * 2
        p1, e1 = mul(ca, cb, N); p2, e2 = mul(sa, sb, N)
        CA, e3 = add(p1, p2, N)
        E_CA = e1 + e2 + e3 + 2 * e_tab * 2
    elif regime == 'small':
        E_red = F(0)
        s = [min(F(1, 2**14), X) + F(1, 2**52)] + xc[1:]
        SAmax = min(F(1), 5 * X)
        CAmax = F(1)
        SBmax = F(1, 2**13)
        sa = canonical(SAmax * ONE_PLUS_U, N); sb = canonical(SBmax * ONE_PLUS_U, N)
        ca = canonical(F(1) * ONE_PLUS_U, N); cb = ca
        e_tab = tab_trunc(N)          # relative to the entry
        p1, e1 = mul(sa, cb, N); p2, e2 = mul(ca, sb, N)
        SA, e3 = add(p1, p2, N)
        E_SA = e1 + e2 + e3 + 2 * e_tab * (SAmax + SBmax)
        p1, e1 = mul(ca, cb, N); p2, e2 = mul(sa, sb, N)
        CA, e3 = add(p1, p2, N)
        E_CA = e1 + e2 + e3 + 2 * e_tab * (1 + SAmax * SBmax)
    else:
        E_red = F(0)
        s = xc
        SAmax, CAmax = F(0), F(1)
        SA = [F(0)] * N
        CA = [F(1)] + [F(0)] * (N - 1)
        E_SA = E_CA = F(0)
    Smax = sum(s, F(0))
    SAv = sum(SA, F(0)); CAv = sum(CA, F(0))

    w, e_w = sqr(s, N)
    Wval = Smax ** 2 + e_w

    def coefb(k, M):
        c = F(1, factorial(k)) * ONE_PLUS_U
        comps = [c * u ** i for i in range(5)]
        return comps[:M], sum(comps[M:], F(0)) + TAB_ERR

    def chain_in(var, Vval, e_var, levels, degs):
        # Horner in var over the coefficients c_k for k in degs (from the top)
        degs = [k for k in degs if k <= (dsin if degs[0] % 2 == 1 else dcos)]
        if not degs:
            return [F(0)] * N, F(0)
        J = len(degs) - 1
        M = levels[J]
        q, e_c = coefb(degs[J], M)
        Eq = e_c
        for j in range(J - 1, -1, -1):
            M = levels[j]
            qM = q + [F(0)] * (M - len(q))
            vM, e_tr = truncate(var, M)
            prod, e_mul = mul(vM, qM, M)
            Eprod = e_mul + e_tr * sum(qM, F(0)) + Vval * Eq
            c, e_c = coefb(degs[j], M)
            q, e_add = add(c, prod, M)
            Eq = e_add + e_c + Eprod
        return q + [F(0)] * (N - len(q)), Eq

    # the chains are evaluated at the computed w (and v): P(w) vs
    # P(s^2) differs by at most |P'| e_w
    dP = poly_deriv(list(range(1, dsin + 1, 2)), Wval)
    dQ = poly_deriv(list(range(2, dcos + 1, 2)), Wval)
    if split == 2:
        P, E_P = chain_in(w, Wval, e_w, lev, list(range(1, dsin + 1, 2)))
        Q, E_Q = chain_in(w, Wval, e_w, lev, list(range(2, dcos + 1, 2)))
        # sin t = SA + (SA w) Q + (CA s) P; cos t = CA + (CA w) Q - (SA s) P
        def combine(Xv, Xb, E_X, Yv, Yb, E_Y):
            m0, e = mul(Xv, w, N); E_m0 = e + E_X * Wval + Xb * e_w
            m1, e = mul(Yv, s, N); E_m1 = e + E_Y * Smax
            l0, e = mul(m0, Q, N); E_l0 = e + E_m0 * sum(Q, F(0)) + sum(m0, F(0)) * (E_Q + dQ * e_w)
            l1, e = mul(m1, P, N); E_l1 = e + E_m1 * sum(P, F(0)) + sum(m1, F(0)) * (E_P + dP * e_w)
            y, e = add3(Xv, l0, l1, N)
            return e + E_X + E_l0 + E_l1
        E_sin = combine(SA, SAv, E_SA, CA, CAv, E_CA)
        E_cos = combine(CA, CAv, E_CA, SA, SAv, E_SA)
    else:
        v, e_v = sqr(w, N)
        Vval = Wval ** 2 + e_v
        P0, E_P0 = chain_in(v, Vval, e_v, lev, list(range(1, dsin + 1, 4)))
        P1, E_P1 = chain_in(v, Vval, e_v, lev, list(range(3, dsin + 1, 4)))
        Q0, E_Q0 = chain_in(v, Vval, e_v, lev, list(range(2, dcos + 1, 4)))
        Q1, E_Q1 = chain_in(v, Vval, e_v, lev, list(range(4, dcos + 1, 4)))
        # the chains at the computed v = w^2 + e_v
        E_P0 += poly_deriv(list(range(1, dsin + 1, 4)), Vval) * e_v
        E_P1 += poly_deriv(list(range(3, dsin + 1, 4)), Vval) * e_v
        E_Q0 += poly_deriv(list(range(2, dcos + 1, 4)), Vval) * e_v
        E_Q1 += poly_deriv(list(range(4, dcos + 1, 4)), Vval) * e_v
        def combine(Xv, Xb, E_X, Yv, Yb, E_Y):
            m0, e = mul(Xv, w, N); E_m0 = e + E_X * Wval + Xb * e_w
            m1, e = mul(m0, w, N); E_m1 = e + E_m0 * Wval + sum(m0, F(0)) * e_w
            m2, e = mul(Yv, s, N); E_m2 = e + E_Y * Smax
            m3, e = mul(m2, w, N); E_m3 = e + E_m2 * Wval + sum(m2, F(0)) * e_w
            ls = []; Es = []
            for (mm, E_mm, ch, E_ch) in [(m0, E_m0, Q0, E_Q0), (m1, E_m1, Q1, E_Q1), (m2, E_m2, P0, E_P0), (m3, E_m3, P1, E_P1)]:
                l, e = mul(mm, ch, N)
                ls.append(l); Es.append(e + E_mm * sum(ch, F(0)) + sum(mm, F(0)) * E_ch)
            y, e = addn([Xv] + ls, N)
            # P = P0 + w P1, Q = Q0 + w Q1 at the computed w
            return e + E_X + sum(Es, F(0)) + sum(m0, F(0)) * dQ * e_w + sum(m2, F(0)) * dP * e_w
        E_sin = combine(SA, SAv, E_SA, CA, CAv, E_CA)
        E_cos = combine(CA, CAv, E_CA, SA, SAv, E_SA)

    # Taylor remainders (|s| <= Smax): sin: s^(dsin+2)/(dsin+2)!, cos: s^(dcos+2)/(dcos+2)!
    Sf = float(Smax)
    R_s = F(Sf ** (dsin + 2) / factorial(dsin + 2)) * F(1001, 1000)
    R_c = F(Sf ** (dcos + 2) / factorial(dcos + 2)) * F(1001, 1000)
    # |sin t - sin t_true| <= E_red, and the remainders enter through
    # CA sin s and SA cos s
    E_sin += CAmax * R_s + SAmax * R_c + E_red + F(1, 2**1000)
    E_cos += SAmax * R_s + CAmax * R_c + E_red + F(1, 2**1000)
    if verbose:
        print("  N = %d split %d regime %s X = %s: sin %s, cos %s (reduction %s)" % (N, split, regime, fr(X), fr(E_sin), fr(E_cos), fr(E_red)))
    return E_sin * (1 + F(1, 2**300)), E_cos * (1 + F(1, 2**300))

def trig_eps(N, split, lev, verbose=False):
    """EPS such that |res - sin x| <= EPS (|x| + min(|x|, 1)) and
    |res - cos x| <= EPS (|x| + 1) for 2^-300 <= |x| <= 2^20 pi/2 - 1,
    checked on a geometric grid (the bounds are monotone in X)."""
    eps = F(0)
    ratio = mpf(2) ** F(1, 4)
    # tiny: 2^-300 <= X < 2^-15 (below 2^-300 the kernels return
    # sin x = x, cos x = 1), bound EPS 2 X for sin and EPS for cos
    X = mpf(2) ** -15
    while X > mpf(2) ** -300:
        Es, Ec = bound_trig(N, split, X, 'tiny', lev)
        eps = max(eps, Es / (2 * X / ratio), Ec / 1)
        X /= ratio
    # small: 2^-15 <= X < 3/4: sin bound EPS 2 X, cos EPS (X + 1)
    X = mpf(2) ** -15 * ratio
    while X < F(3, 4) * ratio:
        Es, Ec = bound_trig(N, split, min(X, F(3, 4)), 'small', lev)
        eps = max(eps, Es / (2 * X / ratio), Ec / (X / ratio + 1))
        X *= ratio
    # general: 3/4 <= X <= XMAX: bounds EPS (X + 1)
    X = F(3, 4) * ratio
    while X < F(1647099) * ratio:
        Es, Ec = bound_trig(N, split, min(X, F(1647099)), 'general', lev)
        eps = max(eps, Es / (X / ratio + 1), Ec / (X / ratio + 1))
        X *= ratio
    if verbose:
        print("N = %d split %d levels %s: TRIG EPS %s" % (N, split, lev, fr(eps)))
    return eps

def bound1_trig(X, regime, verbose=False):
    """The double kernel (N = 1): one table level (A = i/64, |i| <= 50,
    SA = sin A and CA = cos A the nearest doubles), t = fl(h - p1) with
    h = x0 - k P0 and p1 + e1 = k P1 exact, d = t - i/64 exactly
    (|d| <= 1/128), s = fl(d - e1), w = fl(s^2), P and Q by fma chains
    of the rounded coefficients (the constant 1 of P is exact), and
    sin t = fma(SA w, Q, fma(CA s, P, SA)), cos t = fma(CA w, Q,
    fma(-SA s, P, CA)) with the two inner products rounded once each.
    regime 'tiny': X < 1/128 (i = 0: s = x, SA = 0, CA = 1 exactly);
    'small': 1/128 <= X < 3/4 (k = 0, t = x exactly, SA <= 2 X);
    'general': the reduction with |k| <= 2^20.  Returns (E_sin, E_cos)."""
    dsin, dcos = TRIG_DEG[1]
    X = F(X)
    if regime == 'general':
        # |t_true| <= pi/4 + 2^-31 (k is the rounding of fl(x0 c), c the
        # rounded 2/pi: off by one only within 2^-31 of a half-integer),
        # |h - p1 - t_true| = |k| |pi/2 - P0 - P1| <= 2^-67
        H = F(0.7854) + F(1, 2**30)
        E_t = u * H                                  # t = fl(h - p1)
        E1 = u * TRIG_KMAX * abs(PI2[1])             # |e1|
        Sd = F(1, 128) + E1                          # |d - e1|
        S = Sd * ONE_PLUS_U
        E_red = E_t + u * Sd + TRIG_KMAX * PI2_RESID[1]
        SAmax = F(sin(50 / 64)) * F(1001, 1000); CAmax = F(1)
        E_sa = u * SAmax; E_ca = u * CAmax
    elif regime == 'small':
        S = F(1, 128)
        E_red = F(0)
        SAmax = min(F(1), 2 * X); CAmax = F(1)
        E_sa = u * SAmax; E_ca = u * CAmax
    else:
        S = X
        E_red = F(0)
        SAmax = F(0); CAmax = F(1)
        E_sa = E_ca = F(0)
    W = S ** 2 * ONE_PLUS_U
    e_w = u * S ** 2

    def fma_chain(ks):
        # Horner in w from the top: val = fma(w, val, c_k), c_k the
        # rounded 1/k! (exact for k = 1); the error is against the
        # exact polynomial at the computed w
        c = lambda k: F(1, factorial(k))
        val = c(ks[0]) * ONE_PLUS_U
        E = u * c(ks[0])
        for k in ks[1:]:
            newval = (W * val + c(k) * ONE_PLUS_U) * ONE_PLUS_U
            E = u * newval + W * E + (F(0) if k == 1 else u * c(k))
            val = newval
        return val, E
    P, E_P = fma_chain(list(range(dsin, 0, -2)))
    Q, E_Q = fma_chain(list(range(dcos, 1, -2)))
    dP = poly_deriv(list(range(1, dsin + 1, 2)), W)
    dQ = poly_deriv(list(range(2, dcos + 1, 2)), W)
    # the value bounds of the polynomials: |P| <= 1, |Q| <= 1/2 up to
    # the roundings
    def combine(Xb, E_X, Yb, E_Y):
        # fma(X w, Q, fma(Y s, P, X))
        m1 = Yb * S * ONE_PLUS_U; E_m1 = u * Yb * S + E_Y * S
        inner = (m1 * P + Xb) * ONE_PLUS_U
        E_inner = u * inner + E_m1 * P + m1 * (E_P + dP * e_w) + E_X
        m0 = Xb * W * ONE_PLUS_U; E_m0 = u * Xb * W + E_X * W + Xb * e_w
        res = (m0 * Q + inner) * ONE_PLUS_U
        return u * res + E_m0 * Q + m0 * (E_Q + dQ * e_w) + E_inner
    E_sin = combine(SAmax, E_sa, CAmax, E_ca)
    E_cos = combine(CAmax, E_ca, SAmax, E_sa)
    Sf = float(S)
    R_s = F(Sf ** (dsin + 2) / factorial(dsin + 2)) * F(1001, 1000)
    R_c = F(Sf ** (dcos + 2) / factorial(dcos + 2)) * F(1001, 1000)
    E_sin += CAmax * R_s + SAmax * R_c + E_red + F(1, 2**1000)
    E_cos += SAmax * R_s + CAmax * R_c + E_red + F(1, 2**1000)
    if verbose:
        print("  N = 1 regime %s X = %s: sin %s, cos %s (reduction %s)" % (regime, fr(X), fr(E_sin), fr(E_cos), fr(E_red)))
    return E_sin * (1 + F(1, 2**300)), E_cos * (1 + F(1, 2**300))

def trig_eps1(verbose=False):
    """as trig_eps for the double kernel"""
    eps = F(0)
    ratio = mpf(2) ** F(1, 4)
    X = F(1, 128)
    while X > mpf(2) ** -300:
        Es, Ec = bound1_trig(X, 'tiny')
        eps = max(eps, Es / (2 * X / ratio), Ec / 1)
        X /= ratio
    X = F(1, 128) * ratio
    while X < F(3, 4) * ratio:
        Es, Ec = bound1_trig(min(X, F(3, 4)), 'small')
        eps = max(eps, Es / (2 * X / ratio), Ec / (X / ratio + 1))
        X *= ratio
    X = F(3, 4) * ratio
    while X < F(1647099) * ratio:
        Es, Ec = bound1_trig(min(X, F(1647099)), 'general')
        eps = max(eps, Es / (X / ratio + 1), Ec / (X / ratio + 1))
        X *= ratio
    if verbose:
        print("N = 1 (double kernel): TRIG EPS %s" % fr(eps))
    return eps

# ---- the division kernel ----

def div(x, y, Ymin, N):
    """The long division _div_level(q, x, y, N, N): x, y component
    bounds, Ymin a lower bound on |y_0| (so |y| >= Ymin (1 - 2u)).
    Digit k: q_k = fl(r_0 yinv) (|q_k y_0 - r_0| <= (2u + u^2) |r_0|;
    the last digit is fl(r_0 / y_0), within u), the remainder r <-
    (d_0 - e) + (r_1, ...) - q_k (y_1, ..., y_{M-1}) at level M - 1 =
    N - k - 1 (X = TwoSum(d_0, -e), a mul_d and an add3 at that
    level; the components y_M, ... are not subtracted at all).  x/y -
    sum q_k = (r^(N) + dropped) / y with r^(N) the exact final
    remainder.  Returns (digit bounds, |q - x/y|)."""
    r = list(x)
    Qs = []
    err = F(0)
    for k in range(N):
        M = N - k
        delta = u if M == 1 else (2 * u + u * u)
        Qk = r[0] / Ymin * (1 + delta)
        Qs.append(Qk)
        if k + 1 < N:
            X0 = delta * r[0] * ONE_PLUS_U
            X1 = u * X0
            if M - 1 == 1:
                X = [X0]
                err += X1
            else:
                X = [X0, X1] + [F(0)] * (M - 3)
            Y = r[1:M]
            Z, eZ = mul_d(y[1:M], Qk, M - 1)
            err += eZ + Qk * sum(y[M:], F(0))
            r, eA = add3(X, Y, Z, M - 1)
            err += eA
    err += u * r[0] + Qs[-1] * sum(y[1:], F(0))
    return Qs, err / (Ymin * (1 - 2 * u))

def coef_bounds(val, M):
    """the nearest 5-term expansion of |val|, truncated to M terms"""
    c = val * ONE_PLUS_U
    comps = [c * u ** i for i in range(5)]
    return comps[:M], sum(comps[M:], F(0)) + TAB_ERR

def chain_general(var, Vval, levels, coefs):
    """Horner from the top in var (component bounds, value bound
    Vval) over the coefficient magnitudes coefs[j], j = 0..J
    (nearest expansions truncated to levels[j] components; a zero
    coefficient is exact); the error is against the exact polynomial
    at the computed var.  Returns (bounds, value bound, error)."""
    J = len(coefs) - 1
    while J >= 0 and coefs[J] == 0:
        J -= 1
    if J < 0:
        return [F(0)] * len(var), F(0), F(0)
    M = levels[J]
    q, e_c = coef_bounds(coefs[J], M)
    Eq = e_c
    for j in range(J - 1, -1, -1):
        M = levels[j]
        qM = q + [F(0)] * (M - len(q))
        vM, e_tr = truncate(var, M)
        prod, e_mul = mul(vM, qM, M)
        Eprod = e_mul + e_tr * sum(qM, F(0)) + Vval * Eq
        c, e_c = coef_bounds(coefs[j], M) if coefs[j] != 0 else ([F(0)] * M, F(0))
        q, e_add = add(c, prod, M)
        Eq = e_add + e_c + Eprod
    return q + [F(0)] * (len(var) - len(q)), sum(q, F(0)), Eq

def deriv_bound(coefs, Z):
    """|d/dz sum_m coefs[m] z^m| <= this for |z| <= Z"""
    return sum((m * coefs[m] * Z ** (m - 1) for m in range(1, len(coefs))), F(0)) * ONE_PLUS_U

# ---- the logarithm ----

import struct
LOG_M0_BITS = 0x3fe6a00000000000
def from_bits(b):
    return struct.unpack("<d", struct.pack("<Q", b))[0]
def log_interval(ix):
    """[a, b), m_mid and the rounded reciprocal r of the level-1 entry"""
    a = F(from_bits(LOG_M0_BITS + ix * 2**45))
    b = F(from_bits(LOG_M0_BITS + (ix + 1) * 2**45))
    if ix in (74, 75):
        return a, b, F(1), F(1)
    mid = F(from_bits(LOG_M0_BITS + ix * 2**45 + 2**44))
    r = F(float(1 / mid))
    return a, b, mid, r
LOG_DEG = {1: 8, 2: 8, 3: 12, 4: 15}
LOG_LEVELS = {2: [2, 2, 1, 1], 3: [3, 3, 2, 2, 1, 1], 4: [4, 4, 3, 3, 2, 2, 1, 1]}
LOG_LEVELS4 = {2: [2, 1], 3: [3, 2, 1], 4: [4, 3, 2, 1]}

def bound_log(N, split, lev, regime, param, verbose=False, lp=False):
    """The absolute error of dN_log.  regime 'one': e = 0, ix in {74,
    75} (r = 1, z = m - 1 exactly), param = Z an upper bound on |z|
    (and the bound is for |z| in [Z / ratio, Z]); 'far': e = 0, param =
    ix (r != 1); 'exp': param = E = |e| >= 1 (the worst case over the
    level-1 entries).  The reduction: z2 = m (r r') - 1 = (h0 - 1) +
    l0 + tail Rh + m Rl with r r' = Rh + Rl and m0 Rh = h0 + l0 exact,
    the two products with N components, the four terms summed by the
    addn kernel (exact where r = r' = 1; Rl = 0 where r = 1).
    Returns (E_total, R_min), R_min a lower bound on |log x|."""
    deg = LOG_DEG[N]
    Mmax = F(1.4140625)
    m = canonical(Mmax, N)
    if lp:
        # log1p (|x0| <= 1/2): the argument (s, t, x1, ..., x_{N-1})
        # scaled by at most 2, with t <= u s and x_k <= u^k |x0|: the
        # tail (t, x1, ...) and mfull = (s, t, x1, ..., x_{N-2}); |x0|
        # is at most Z + u in the regime 'one' (e = 0), max |m - 1| in
        # 'far' (e = 0), and 1/2 in 'exp' (then |e| = 1)
        if regime == 'one':
            Xb = F(param) + u
        elif regime == 'far':
            a, b, _, _ = log_interval(param)
            Xb = max(abs(a - 1), abs(b - 1))
        else:
            assert param == 1
            Xb = F(1, 2)
        tail = [u * Mmax] + [u ** k * Xb * 2 for k in range(1, N)]
        m = [Mmax] + tail[:N - 1]
    else:
        tail = m[1:] + [F(0)]
    Rhmax = F(1) / F(0.70703125) * (1 + F(1, 2**7)) * ONE_PLUS_U
    Z2true = F(1, 2**14) * (1 + F(1, 2**7)) + 3 * u     # |m r r' - 1|
    eL, e_eL, c0 = [F(0)] * N, F(0), F(0)
    if regime == 'one':
        Z = F(param)
        Rmin = Z / (mpf(2) ** F(1, 4)) * (1 - Z / 2)
        Lb, e_L = [F(0)] * N, F(0)
        rp1 = min((Z + F(1, 2**14)) * ONE_PLUS_U + u, F(1, 2**7) + u)   # |r' - 1|
        if Z <= F(1, 2**14):
            # r = r' = 1: z2 = z = m - 1 exactly
            Z2 = Z; z2 = canonical(Z, N); E_z2 = F(0)
            Lp, e_Lp = [F(0)] * N, F(0)
        else:
            # r = 1, Rh = r', Rl = 0: z2 = (h0 - 1) + l0 + tail r'
            Rh = 1 + rp1
            H = Z2true + 2 * u * Mmax * Rh
            c = [H * ONE_PLUS_U + u] + [F(0)] * (N - 1)
            D = [u * (1 + H)] + [F(0)] * (N - 1)
            P, e_P = mul_d(tail, Rh, N)
            z2, e_a = add3(c, D, P, N)
            E_z2 = e_a + e_P
            Z2 = Z2true + E_z2
            z2 = canonical(Z2, N)
            Lp, e_Lp = coef_bounds(rp1 * (1 + rp1), N)
    else:
        if regime == 'far':
            ix = param
            a, b, mid, r = log_interval(ix)
            Lval = abs(mp.log(r))
            Rmin = min(abs(mp.log(a)), abs(mp.log(b)))
        else:
            E = F(param)
            Lval = abs(mp.log(F(0.70703125)))
            Rmin = E * mp.log(2) - abs(mp.log(F(0.70703125)))
            # e L0 exact, e (L1 + ... + LN) by mul_d, plus the residual
            c0 = E * abs(L[0])
            Lt = [abs(L[i]) for i in range(1, N + 1)]
            eL, e_eL = mul_d(Lt, E, N)
            e_eL += E * LN2_RESID[N]
        Lb, e_L = coef_bounds(Lval, N)
        rp1 = F(1, 2**7) + u
        Lp, e_Lp = coef_bounds(rp1 * (1 + rp1), N)
        Rlmax = u * Rhmax
        H = Z2true + Rlmax * Mmax + 2 * u * Mmax * Rhmax
        c = [H * ONE_PLUS_U + u] + [F(0)] * (N - 1)
        D = [u * (1 + H)] + [F(0)] * (N - 1)
        P, e_P = mul_d(tail, Rhmax, N)
        Q, e_Q = mul_d(m, Rlmax, N)
        z2, e_a = addn([c, D, P, Q], N)
        E_z2 = e_a + e_P + e_Q
        Z2 = Z2true + E_z2
        z2 = canonical(Z2, N)
    # D = c0 + eL + L + L' + z2
    c0v = [c0] + [F(0)] * (N - 1)
    Dv, e_D = addn([c0v, eL, Lb, Lp, z2], N)
    E_D = e_D + e_eL + e_L + e_Lp + E_z2
    Dval = sum(Dv, F(0))
    # the polynomial
    w, e_w = sqr(z2, N)
    Wval = Z2 ** 2 + e_w
    cE = [F(1, k) for k in range(2, deg + 1, 2)]     # E(w): z^2, z^4, ...
    cO = [F(1, k) for k in range(3, deg + 1, 2)]     # O(w): z^3, z^5, ...
    dE = deriv_bound(cE, Wval); dO = deriv_bound(cO, Wval)
    if split == 2:
        qE, Ev, E_E = chain_general(w, Wval, lev, cE)
        qO, Ov, E_O = chain_general(w, Wval, lev, cO)
        zw, e = mul(z2, w, N); E_zw = e + E_z2 * Wval + Z2 * e_w
        l0, e = mul(w, qE, N); E_l0 = e + e_w * Ev + Wval * (E_E + dE * e_w)
        l1, e = mul(zw, qO, N); E_l1 = e + E_zw * Ov + sum(zw, F(0)) * (E_O + dO * e_w)
        y, e = add3(Dv, l0, l1, N)
        E_total = e + E_D + E_l0 + E_l1
    else:
        v, e_v = sqr(w, N)
        Vval = Wval ** 2 + e_v
        chains = []
        for cs in (cE[0::2], cE[1::2], cO[0::2], cO[1::2]):
            q, val, err = chain_general(v, Vval, lev, cs)
            err += deriv_bound(cs, Vval) * e_v
            chains.append((q, val, err))
        zw, e = mul(z2, w, N); E_zw = e + E_z2 * Wval + Z2 * e_w
        zv, e = mul(z2, v, N); E_zv = e + E_z2 * Vval + Z2 * e_v
        mults = [(w, e_w), (v, e_v), (zw, E_zw), (zv, E_zv)]
        ls = []; Es = F(0)
        for (mm, E_mm), (q, val, err) in zip(mults, chains):
            l, e = mul(mm, q, N)
            ls.append(l); Es += e + E_mm * val + sum(mm, F(0)) * err
        # E = E0 + w E1, O = O0 + w O1 at the computed w
        Es += Wval * dE * e_w + sum(zw, F(0)) * dO * e_w
        y, e = addn([Dv] + ls, N)
        E_total = e + E_D + Es
    # the polynomial remainder beyond degree deg (|z2| <= Z2)
    R = Z2 ** (deg + 1) / (deg + 1) / (1 - Z2) * F(1001, 1000)
    E_total += R + F(1, 2**1000)
    if verbose:
        print("  log N = %d split %d regime %s param %s: error %s, |log x| >= %s, rel %s" % (N, split, regime, param, fr(E_total), fr(Rmin), fr(E_total / Rmin)))
    return E_total * (1 + F(1, 2**300)), Rmin

def log_eps(N, split, lev, verbose=False, lp=False):
    # below |z2| = 2^-250 the kernel returns log1p(z2) = z2: a relative
    # error of at most 2^-251
    eps = mpf(2) ** -251
    ratio = mpf(2) ** F(1, 4)
    Z = F(1, 2**7)
    while Z > mpf(2) ** -250:
        E, R = bound_log(N, split, lev, 'one', Z, lp=lp)
        eps = max(eps, E / R)
        Z /= ratio
    for ix in list(range(0, 74)) + list(range(76, 128)):
        E, R = bound_log(N, split, lev, 'far', ix, lp=lp)
        eps = max(eps, E / R)
    E_ = 1
    while E_ <= (1 if lp else 2048):
        E, R = bound_log(N, split, lev, 'exp', E_, lp=lp)
        eps = max(eps, E / R)
        E_ *= 2
    if lp:
        # |x0| > 1/2: 1 + x rounded to N components by the add kernel
        # (a relative error e/|1 + x|, absolute in the log) and the
        # plain kernel; |log(1 + x)| >= log(3/2)
        e_log = log_eps(N, split, lev)
        worst = F(0)
        X = F(1, 2)
        while X < mpf(2) ** 1024:
            one = [F(1)] + [F(0)] * (N - 1)
            y, e = add(one, canonical(X, N), N)
            worst = max(worst, e / (X - 1) if X > 1 else e / F(1, 2))   # against |1 + x| >= X - 1 (or 1/2)
            X *= 2
        eps = max(eps, e_log * (1 + worst) + worst / mp.log(F(3, 2)))
    if verbose:
        print("N = %d split %d levels %s: %s EPS %s" % (N, split, lev, "LOG1P" if lp else "LOG", fr(eps)))
    return eps

def bound1_log(regime, param, verbose=False, lp=False):
    """The double kernel: z = fma(m, r, -1) + t r (a relative error u,
    plus the roundings of t r and of the sum for log1p, where t <= u m
    is the low part of the argument 1 + x; exact where r = 1 and t =
    0), the polynomial of degree 8 by fma with the rounded
    coefficients, d = ((e L0 + L) + e L1) + z rounded at each step
    and res = fma(w, P, d).  regime as for bound_log."""
    deg = 8
    Mmax = F(1.4140625)
    tr = (u * Mmax * F(1) / F(0.70703125) * (1 + u)) if lp else F(0)   # |t r|
    if regime == 'one':
        Z = F(param); Rmin = Z / (mpf(2) ** F(1, 4)) * (1 - Z / 2)
        # z = (m - 1) + t exactly summed then rounded: a relative error u
        E_z = u * Z * ONE_PLUS_U if lp else F(0); Lval = F(0); c0 = F(0); E = F(0)
    elif regime == 'far':
        a, b, mid, r = log_interval(param)
        Z = (max(abs(a * r - 1), abs(b * r - 1)) + tr) * ONE_PLUS_U
        Rmin = min(abs(mp.log(a)), abs(mp.log(b)))
        E_z = 2 * u * Z + u * tr; Lval = abs(mp.log(r)); c0 = F(0); E = F(0)
    else:
        E = F(param); Z = (F(1, 2**7) + tr) * ONE_PLUS_U
        Lval = abs(mp.log(F(0.70703125))); Rmin = E * mp.log(2) - Lval
        E_z = 2 * u * Z + u * tr; c0 = E * abs(L[0])
    W = Z ** 2 * ONE_PLUS_U; e_w = u * Z ** 2
    # P = fma(w, P, fma(z, c_odd, c_even)) from the top: c8; (c6, c7); (c4, c5); (c2, c3)
    cs = [F(1, k) for k in range(2, 9)]   # |c_2| .. |c_8|
    val = cs[6] * ONE_PLUS_U; Eq = u * cs[6]
    for (ce, co) in ((cs[4], cs[5]), (cs[2], cs[3]), (cs[0], cs[1])):
        inner = (Z * co + ce) * ONE_PLUS_U
        E_inner = u * inner + u * co * Z + E_z * co + u * ce
        newval = (W * val + inner) * ONE_PLUS_U
        Eq = u * newval + W * Eq + val * e_w + E_inner
        val = newval
    # P(w, z) = sum_k (c_{2k+2} + c_{2k+3} z) w^k at the computed w
    Eq += (deriv_bound(cs[0::2], W) + Z * deriv_bound(cs[1::2], W)) * e_w
    # d = ((e L0 + L) + e L1) + z
    d1 = (c0 + Lval) * ONE_PLUS_U; E_d1 = u * d1 + u * Lval
    d2 = (d1 + E * abs(L[1])) * ONE_PLUS_U; E_d2 = u * d2 + E_d1 + E * LN2_RESID[1]
    d3 = (d2 + Z) * ONE_PLUS_U; E_d3 = u * d3 + E_d2 + E_z
    res = (W * val + d3) * ONE_PLUS_U
    E_total = u * res + W * Eq + val * e_w + E_d3
    R = Z ** (deg + 1) / (deg + 1) / (1 - Z) * F(1001, 1000)
    E_total += R + F(1, 2**1000)
    if verbose:
        print("  log N = 1 regime %s param %s: error %s, rel %s" % (regime, param, fr(E_total), fr(E_total / Rmin)))
    return E_total * (1 + F(1, 2**300)), Rmin

def log_eps1(verbose=False, lp=False):
    # (below |z| = 2^-500 the double kernel returns z: w = 0)
    eps = F(0)
    ratio = mpf(2) ** F(1, 4)
    Z = F(1, 2**7)
    while Z > mpf(2) ** -500:
        E, R = bound1_log('one', Z, lp=lp)
        eps = max(eps, E / R)
        Z /= ratio
    for ix in list(range(0, 74)) + list(range(76, 128)):
        E, R = bound1_log('far', ix, lp=lp)
        eps = max(eps, E / R)
    E_ = 1
    while E_ <= 2048:
        E, R = bound1_log('exp', E_, lp=lp)
        eps = max(eps, E / R)
        E_ *= 2
    if verbose:
        print("N = 1 (double kernel): %s EPS %s" % ("LOG1P" if lp else "LOG", fr(eps)))
    return eps

# ---- the arctangent ----

ATAN_K = {1: 2, 2: 6, 3: 10, 4: 13}
ATAN_LEVELS = {2: [2, 2, 1, 1], 3: [3, 3, 2, 2, 1, 1], 4: [4, 4, 3, 3, 2, 2, 1]}
ATAN_LEVELS8 = {2: [2, 1], 3: [3, 2, 1], 4: [4, 3, 2, 1]}

def bound_atan(N, split, lev, regime, param, verbose=False):
    """The absolute error of dN_atan for x >= 0.  regime 'zero': i =
    0, param = X <= 1/128 (the bound is for x in [X / ratio, X]);
    'small': x <= 1, param = i >= 1; 'big': x > 1, param = (Xmin,
    Xmax) the range of x (i = round(64 / x) must be constant on it);
    'huge': x >= 2^900, param = X.  Returns (E_total, R_min)."""
    K = ATAN_K[N]
    if regime == 'huge':
        X = F(param)
        # pi/2 - fl(1/x0): 1/x vs 1/x0 (the tail of x), the rounding,
        # the neglected 1/(3 x^3), and the add with the 5-term pi/2
        t = [u / X * 2 + F(1, X)] + [F(0)] * (N - 1)   # |t| and its error below
        E_t = u / X * 2 + u * F(1, X) + F(1, 3 * X ** 3)
        Pb, e_P = coef_bounds(mp.pi / 2, N)
        y, e = add(Pb, t, N)
        E_total = e + e_P + E_t + F(1, 2**1000)
        Rmin = mp.pi / 4
        return E_total * (1 + F(1, 2**300)), Rmin
    if regime == 'zero':
        X = F(param)
        z = canonical(X, N); Zc = X; E_z = F(0)
        Cb, e_C = [F(0)] * N, F(0)
        Rmin = X / (mpf(2) ** F(1, 4)) * (1 - X ** 2 / 3)
    elif regime == 'small':
        i = param
        c = F(i, 64)
        xmax = min(F(1), c + F(1, 128)); xmin = c - F(1, 128)
        x = canonical(xmax, N)
        # n = x - c: the head exact, the tail that of x
        n = [F(1, 128)] + x[1:]
        n = canonical(sum(n, F(0)), N)
        t, e_t = mul_d(x, c, N)
        one = [F(1)] + [F(0)] * (N - 1)
        d, e_d = add(one, t, N)
        Ymin = (1 + c * xmin) * (1 - 2 * u)
        Qs, e_div = div(n, d, Ymin, N)
        Zc = sum(Qs, F(0))
        E_z = e_div + Zc * (e_t + e_d) / Ymin
        z = canonical(Zc, N)
        Cb, e_C = coef_bounds(mp.atan(c), N)
        Rmin = mp.atan(xmin)
    else:
        Xmin, Xmax = F(param[0]), F(param[1])
        i = int(mp.nint(64 / Xmin))
        assert i == int(mp.nint(64 / Xmax)), (Xmin, Xmax)
        c = F(i, 64)
        x = canonical(Xmax, N)
        t, e_t = mul_d(x, c, N)
        # n = 1 - c x: the head exact (c x in [1/2, 2] when i >= 1; c = 0 otherwise)
        Nmax = max(abs(1 - c * Xmin), abs(1 - c * Xmax)) + 2 * u * (1 + c * Xmax)
        n = [Nmax * ONE_PLUS_U] + t[1:]
        n = canonical(sum(n, F(0)), N)
        cc = [c] + [F(0)] * (N - 1)
        d, e_d = add(x, cc, N)
        Ymin = (Xmin + c) * (1 - 2 * u)
        Qs, e_div = div(n, d, Ymin, N)
        Zc = sum(Qs, F(0))
        E_z = e_div + Zc * e_d / Ymin + e_t / Ymin
        z = canonical(Zc, N)
        Pb, e_P = coef_bounds(mp.pi / 2, N)
        Ab, e_A = coef_bounds(mp.atan(c), N)
        Cb, e = add(Pb, Ab, N)
        e_C = e + e_P + e_A
        Rmin = mp.atan(Xmin)
    # D = C + z
    Dv, e = add(Cb, z, N)
    E_D = e + e_C + E_z
    w, e_w = sqr(z, N)
    Wval = Zc ** 2 + e_w
    cE = [F(1, 2 * k + 3) for k in range(0, K + 1)]   # |E| coefficients of w^k
    dE = deriv_bound(cE, Wval)
    p, e = mul(z, w, N); E_p = e + E_z * Wval + Zc * e_w
    pval = sum(p, F(0))
    if split == 2:
        v, e_v = sqr(w, N)
        Vval = Wval ** 2 + e_v
        q0, v0, E_0 = chain_general(v, Vval, lev, cE[0::2])
        q1, v1, E_1 = chain_general(v, Vval, lev, cE[1::2])
        E_0 += deriv_bound(cE[0::2], Vval) * e_v
        E_1 += deriv_bound(cE[1::2], Vval) * e_v
        l0, e = mul(p, q0, N); E_l0 = e + E_p * v0 + pval * E_0
        pw, e = mul(p, w, N); E_pw = e + E_p * Wval + pval * e_w
        l1, e = mul(pw, q1, N); E_l1 = e + E_pw * v1 + sum(pw, F(0)) * E_1
        y, e = add3(Dv, l0, l1, N)
        E_total = e + E_D + E_l0 + E_l1 + pval * dE * e_w
    else:
        w2, e_w2 = sqr(w, N)
        W2val = Wval ** 2 + e_w2
        uu, e_u = sqr(w2, N)
        Uval = W2val ** 2 + e_u
        w3, e = mul(w, w2, N); E_w3 = e + e_w * W2val + Wval * e_w2
        # the multipliers p (1, w, w2, w3): the first exact
        mults = [(p, E_p)]
        for (mm, E_mm) in ((w, e_w), (w2, e_w2), (w3, E_w3)):
            l, e = mul(p, mm, N)
            mults.append((l, e + E_p * sum(mm, F(0)) + pval * E_mm))
        ls = []; Es = F(0)
        for lane in range(4):
            cs = cE[lane::4]
            q, val, err = chain_general(uu, Uval, lev, cs)
            err += deriv_bound(cs, Uval) * e_u
            mm, E_mm = mults[lane]
            l, e = mul(mm, q, N)
            ls.append(l); Es += e + E_mm * val + sum(mm, F(0)) * err
        # E(w) = sum w^k E_k(w^4) at the computed w
        Es += pval * dE * e_w
        y, e = addn([Dv] + ls, N)
        E_total = e + E_D + Es
    R = Zc ** (2 * K + 5) / (2 * K + 5) / (1 - Zc ** 2) * F(1001, 1000)
    E_total += R + F(1, 2**1000)
    if verbose:
        print("  atan N = %d split %d regime %s param %s: error %s, |atan x| >= %s, rel %s" % (N, split, regime, param, fr(E_total), fr(Rmin), fr(E_total / Rmin)))
    return E_total * (1 + F(1, 2**300)), Rmin

def atan_eps(N, split, lev, verbose=False):
    eps = F(0)
    ratio = mpf(2) ** F(1, 4)
    X = F(1, 128)
    while X > mpf(2) ** -200:
        E, R = bound_atan(N, split, lev, 'zero', X)
        eps = max(eps, E / R)
        X /= ratio
    for i in range(1, 65):
        E, R = bound_atan(N, split, lev, 'small', i)
        eps = max(eps, E / R)
    # x > 1: a grid fine enough for i = round(64/x) to be constant
    # (the breakpoints are 64/(i + 1/2)); beyond 128, i = 0
    pts = sorted(set([F(1)] + [F(128, 2 * i + 1) for i in range(0, 64)] + [F(128)]))
    for a, b in zip(pts[:-1], pts[1:]):
        # subdivide by the ratio
        lo = a
        while lo < b:
            hi = min(lo * ratio, b)
            # avoid the breakpoint itself on the wrong side
            E, R = bound_atan(N, split, lev, 'big', (lo * (1 + F(1, 2**60)), hi * (1 - F(1, 2**60))))
            eps = max(eps, E / R)
            lo = hi
    X = F(128)
    while X < mpf(2) ** 900:
        E, R = bound_atan(N, split, lev, 'big', (X, X * 2))
        eps = max(eps, E / R)
        X *= 2
    E, R = bound_atan(N, split, lev, 'huge', mpf(2) ** 900)
    eps = max(eps, E / R)
    if verbose:
        print("N = %d split %d levels %s: ATAN EPS %s" % (N, split, lev, fr(eps)))
    return eps

def bound1_atan(regime, param, verbose=False):
    """The double kernel: z = (x - c) / fl(1 + x c) or fl(1 - c x) /
    (x + c) (two or three roundings), E = fma chain of three
    coefficients, res = C + fma(z w, E, z)."""
    K = 2
    if regime == 'zero':
        X = F(param); Z = X; E_z = F(0); C = F(0); E_C = F(0)
        Rmin = X / (mpf(2) ** F(1, 4)) * (1 - X ** 2 / 3)
    elif regime == 'small':
        i = param; c = F(i, 64)
        xmin = c - F(1, 128)
        Z = F(1, 128) / (1 + c * xmin) * (1 + 3 * u)
        E_z = 3 * u * Z
        C = mp.atan(c); E_C = u * C
        Rmin = mp.atan(xmin)
    elif regime == 'big':
        Xmin, Xmax = F(param[0]), F(param[1])
        i = int(mp.nint(64 / Xmin)); c = F(i, 64)
        Nmax = max(abs(1 - c * Xmin), abs(1 - c * Xmax))
        Z = Nmax / (Xmin + c) * (1 + 4 * u)
        E_z = 4 * u * Z + u * (1 + c * Xmax) / (Xmin + c)    # the rounding of c x is absolute
        C = mp.pi / 2 - mp.atan(c); E_C = u * mp.pi / 2 + u * mp.atan(c) + u * C
        Rmin = mp.atan(Xmin)
    else:
        X = F(param)
        # pi/2 - fl(1/x)
        E_total = u * mp.pi / 2 + u / X + F(1, 3 * X ** 3) + u * mp.pi / 2
        return E_total * (1 + F(1, 2**300)), mp.pi / 4
    W = Z ** 2 * ONE_PLUS_U; e_w = u * Z ** 2
    cs = [F(1, 2 * k + 3) for k in range(0, K + 1)]
    val = cs[K] * ONE_PLUS_U; Eq = u * cs[K]
    for k in range(K - 1, -1, -1):
        newval = (W * val + cs[k]) * ONE_PLUS_U
        Eq = u * newval + W * Eq + val * e_w + u * cs[k]
        val = newval
    Eq += deriv_bound(cs, W) * e_w
    p = Z * W * ONE_PLUS_U; E_p = u * Z * W + E_z * W + Z * e_w
    inner = (p * val + Z) * ONE_PLUS_U
    E_inner = u * inner + E_p * val + p * Eq + E_z
    res = (C + inner) * ONE_PLUS_U
    E_total = u * res + E_C + E_inner
    R = Z ** (2 * K + 5) / (2 * K + 5) / (1 - Z ** 2) * F(1001, 1000)
    E_total += R + F(1, 2**1000)
    if verbose:
        print("  atan N = 1 regime %s param %s: error %s, rel %s" % (regime, param, fr(E_total), fr(E_total / Rmin)))
    return E_total * (1 + F(1, 2**300)), Rmin

def atan_eps1(verbose=False):
    eps = F(0)
    ratio = mpf(2) ** F(1, 4)
    X = F(1, 128)
    while X > mpf(2) ** -200:
        E, R = bound1_atan('zero', X)
        eps = max(eps, E / R)
        X /= ratio
    for i in range(1, 65):
        E, R = bound1_atan('small', i)
        eps = max(eps, E / R)
    pts = sorted(set([F(1)] + [F(128, 2 * i + 1) for i in range(0, 64)] + [F(128)]))
    for a, b in zip(pts[:-1], pts[1:]):
        lo = a
        while lo < b:
            hi = min(lo * ratio, b)
            E, R = bound1_atan('big', (lo * (1 + F(1, 2**60)), hi * (1 - F(1, 2**60))))
            eps = max(eps, E / R)
            lo = hi
    X = F(128)
    while X < mpf(2) ** 900:
        E, R = bound1_atan('big', (X, X * 2))
        eps = max(eps, E / R)
        X *= 2
    E, R = bound1_atan('huge', mpf(2) ** 900)
    eps = max(eps, E / R)
    if verbose:
        print("N = 1 (double kernel): ATAN EPS %s" % fr(eps))
    return eps

def eps_const(name, N, te):
    e = int(mp.floor(mp.log(te, 2))); m = int(mp.ceil(te / mpf(2) ** (e - 7)))
    const = float(m * mpf(2) ** (e - 7))
    assert const >= te
    print("  DFLOAT_%s_EPS_%d %s\n" % (name, N, const.hex()))

TRIG_LEVELS = {1: [1, 1, 1, 1], 2: [2, 2, 1, 1], 3: [3, 3, 2, 1, 1, 1], 4: [4, 4, 3, 2, 2, 1, 1]}
TRIG_LEVELS4 = {1: [1, 1], 2: [2, 1], 3: [3, 2, 1], 4: [4, 3, 2, 1]}

def trig_search(N, split):
    import itertools
    dsin, dcos = TRIG_DEG[N]
    J = (dsin - 1) // 2 if split == 2 else ((dsin - 1) // 2 + 2) // 2 - 1
    full = trig_eps(N, split, [N] * (J + 1))
    best = None
    for tail in itertools.combinations_with_replacement(range(1, N + 1), J):
        lev = [N] + sorted(tail, reverse=True)
        b = trig_eps(N, split, lev)
        if b <= 2 * full:
            cost = sum(m * m for m in lev)
            if best is None or cost < best[0]:
                best = (cost, lev, b)
    print("trig N = %d split %d: full-precision %s; cheapest within 2x: %s, %s" % (N, split, fr(full), best[1], fr(best[2])))
    return best

def search(N, split=2):
    """the cheapest schedule (fewest components in total) whose bound is
    within a factor 2 of the full-precision one"""
    deg = {1: 3, 2: 6, 3: 8, 4: 11}[N]
    JE = (deg - 1) // 2
    J = JE if split == 2 else (JE + 2) // 2 - 1
    full = bound(N, [N] * (J + 1), split)
    best = None
    import itertools
    for tail in itertools.combinations_with_replacement(range(1, N + 1), J):
        lev = [N] + sorted(tail, reverse=True)
        b = bound(N, lev, split)
        if b <= 2 * full:
            cost = sum(m * m for m in lev)
            if best is None or cost < best[0]:
                best = (cost, lev, b)
    print("N = %d, split %d: full-precision bound %s; cheapest within 2x: %s, bound %s" % (N, split, fr(full), best[1], fr(best[2])))
    return best

# the schedules used in template.inc (from --search: the cheapest
# within a factor 2 of the full-precision bound)
LEVELS = {
    1: [1, 1],
    2: [2, 1, 1],
    3: [3, 2, 2, 1],
    4: [4, 3, 3, 2, 1, 1],
}
LEVELS4 = {
    1: [1],
    2: [2, 1],
    3: [3, 2],
    4: [4, 3, 1],
}

if __name__ == "__main__":
    if "--search" in sys.argv:
        for N in range(2, 5):
            search(N, 2)
            search(N, 4)
    elif "--trig-search" in sys.argv:
        for N in range(2, 5):
            trig_search(N, 2)
            trig_search(N, 4)
    else:
        for N in range(1, 5):
            if N == 1:
                rel = bound1(verbose=True)
            else:
                rel2 = bound(N, LEVELS[N], 2, verbose=True)
                rel4 = bound(N, LEVELS4[N], 4, verbose=True)
                rel = max(rel2, rel4)
            if "--expm1" in sys.argv:
                if N == 1:
                    te = expm1_eps1(verbose=True)
                else:
                    # the two-way (non-SIMD) path runs the chains untapered for expm1
                    t2 = expm1_eps(N, 2, [N] * len(LEVELS[N]), verbose=True)
                    t4 = expm1_eps(N, 4, LEVELS4[N], verbose=True)
                    te = max(t2, t4)
                eps_const("EXPM1", N, te)
                continue
            if "--log1p" in sys.argv:
                if N == 1:
                    te = log_eps1(verbose=True, lp=True)
                else:
                    t2 = log_eps(N, 2, LOG_LEVELS[N], verbose=True, lp=True)
                    t4 = log_eps(N, 4, LOG_LEVELS4[N], verbose=True, lp=True)
                    te = max(t2, t4)
                eps_const("LOG1P", N, te)
                continue
            if "--log" in sys.argv or "--atan" in sys.argv:
                if "--log" in sys.argv:
                    if N == 1:
                        te = log_eps1(verbose=True)
                    else:
                        t2 = log_eps(N, 2, LOG_LEVELS[N], verbose=True)
                        t4 = log_eps(N, 4, LOG_LEVELS4[N], verbose=True)
                        te = max(t2, t4)
                    eps_const("LOG", N, te)
                if "--atan" in sys.argv:
                    if N == 1:
                        te = atan_eps1(verbose=True)
                    else:
                        t2 = atan_eps(N, 2, ATAN_LEVELS[N], verbose=True)
                        t8 = atan_eps(N, 8, ATAN_LEVELS8[N], verbose=True)
                        te = max(t2, t8)
                    eps_const("ATAN", N, te)
                continue
            if "--trig" in sys.argv:
                if N == 1:
                    te = trig_eps1(verbose=True)
                else:
                    t2 = trig_eps(N, 2, TRIG_LEVELS[N], verbose=True)
                    t4 = trig_eps(N, 4, TRIG_LEVELS4[N], verbose=True)
                    te = max(t2, t4)
                e = int(mp.floor(mp.log(te, 2))); m = int(mp.ceil(te / mpf(2) ** (e - 7)))
                const = float(m * mpf(2) ** (e - 7))
                assert const >= te
                print("  DFLOAT_TRIG_EPS_%d %s\n" % (N, const.hex()))
                continue
            # the bound as a double with 8 significant bits, rounded up
            e = int(mp.floor(mp.log(rel, 2)))
            m = int(mp.ceil(rel / mpf(2) ** (e - 7)))
            const = float(m * mpf(2) ** (e - 7))
            assert const >= rel
            print("  DFLOAT_EXP_EPS_%d %s\n" % (N, const.hex()))
