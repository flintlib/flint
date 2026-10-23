#!/usr/bin/env python3
"""Emit src/mp_real/dense_tables.c and src/mp_real/dense_tables.h: the
tables of the dense-table per-size kernels (the _mp_real_*_opt_<n>
files that dev/tune_mp_real.py emits for n in the dense range, see
dev/gen_mp_real_dense_bodies.py), 64-bit limbs only.

Every table stores its values at W = MP_REAL_DENSE_W fraction limbs,
least significant first, as floors; a kernel at n <= W fraction limbs
reads the top n limbs of an entry (offset W - n), which is the floor at
n limbs.  The first two levels (L_1, G_2, M_1, M_2, A_1, A_2 and the
sin/cos T_1 below) are not emitted here: the kernels read them from the
21-limb tables of dev/gen_mp_real_dense_rt_tables.py, which they share
with the run-time kernels of dense_rt.c.  The reductions consume seven bits per level.

  exp:  C_1[j] = floor((exp(j/2^7) / 2^s - 1) 2^64)      j < 128 (one limb),
                s = 1 for j >= S_1 (exp(j/2^7) >= 2), else 0
        L_1[j] = s log 2 + log(1 + C_1[j] / 2^64)       j < 128
        G_k[j] = log(1 + j 2^-7k)                       j <= J_k, k = 2 .. EXP_LEVELS
        G_top[j] = log(1 + j 2^-64)                     j <= 3 (the top level)

The exp reduction (v in [0, 1) at n limbs) subtracts logarithms of
exact factors at every level: t = v - L_1[j_1] (j_1 the top seven bits;
L_1[j_1] <= j_1/2^7 < L_1[j_1] + 2^-63), then for k = 2, 3, ...:
j_k = floor(t 2^7k) (the bits below 2^-7(k-1) of the top limb, which may
exceed 127 by the overshoot of the previous level), t <- t - G_k[j_k].
As every table value is a floor of log(1 + c) <= c (c = j_k 2^-7k <= t),
t stays nonnegative, and

    t < 2^-7k + max (c - log(1 + c)) + (k - 1) ulps,

which bounds j_{k+1} (J_{k+1} below, computed here with a margin) and,
after the last level L >= 2, the remainder: t < 2^-(7L - 1).  Then
exp(v) = 2^s (1 + C_1[j_1]/2^64) prod (1 + j_k 2^-7k) exp(t), every
factor a single limb.  After all EXP_LEVELS = 9 levels, t < (1 + 2^-8)
2^-63, and the top level takes the whole top limb as its index (j <= 2,
the factor 1 + j 2^-64), leaving t below 2^-64 (unless the next limb is
within a few units of 2^64; then below 2^-63) for the hardcoded series
of the 2^-64 tier.

  log1p: D_1[j] = floor(2^64 j / (2^7 + j)) / 2^64      j < 128 (one limb)
         M_1[j] = -log(1 - D_1[j])                      j < 128
         M_k[j] = -log(1 - D_k(j))                      j <= K_k, k = 2 .. LOG_LEVELS
         with D_k(j) = c - ceil(c^2 2^64) / 2^64, c = j 2^-7k, which
         the kernel forms itself (one limb; c^2 2^64 is exact up to
         k = 4, below 2^-6 beyond)

The log1p reduction (u = v in [0, 1) at n limbs): for k = 1, 2, ...:
j_k = floor(u 2^7k), u <- (1 + u)(1 - D) - 1 = u - D - u D with u D
truncated (so that u stays an upper bound of the exact value, by
under one ulp per level), and acc += M_k[j_k].  Since D <= c / (1 + c)
(for k >= 2, c - c^2 <= c / (1 + c)), u stays nonnegative, and

    u < max_j ((1 + (j + 1) 2^-7k)(1 - D_k(j)) - 1) + k ulps

bounds j_{k+1} (K_{k+1}) and, after level L, u < 2^-(7L - 1).  Then
log(1 + v) = acc + log(1 + u).

  atan: A_1[j] = atan(j/2^7)                            j < 128
        A_k[j] = atan(j 2^-7k)                          j <= N_k, k = 2 .. ATAN_LEVELS

The atan reduction rotates (X, Y) = (1, v) by the exact factors
1 - i c_k, c_k = j_k 2^-7k, which subtract atan(c_k) from the angle:
X <- X + c Y, Y <- Y - c X.  Level 1 takes j_1 = floor(2^7 v) (X = 1);
from then on the index comes from q = Y/X estimated with a one-limb
reciprocal R of X after level 1, the estimate lowered by 2^-12 so that
c_k <= Y/X (X only grows, by under 2^-13.8 relative after level 1):
the angle stays nonnegative and q < 1.042 2^-7k after level k, which
bounds N_k (j <= 133) and leaves q < 2^-(7L - 1) for the series after
the division q = Y/X.

  sin, cos: C_1[j] = floor(tan(j/2^7) 2^64)                 j < 64 (one limb)
            T_1[j] = atan(C_1[j] / 2^64)                    j < 64
            (and the atan tables A_k for k >= 2)

The sin/cos reduction works on w = v/2 < 1/2: t = w - T_1[j_1] (j_1 the
top seven bits of w), then t <- t - A_k[j_k], j_k = floor(t 2^7k), for
k = 2 .. L.  Every subtracted angle is atan(c) <= c <= t for the exact
multiplier c (c_1 = C_1[j_1] / 2^64 <= tan(j_1/2^7)), and t < 2^-7k (1 +
2^-20) after level k, so that j_k <= 128 and the series takes t <
2^-(7L - 1).  Then (cos w, sin w) is (cos t, sin t) rotated by the
factors 1 + i c_k, scaled by sqrt(N), N = prod (1 + c_k^2), and the
double-angle formulas cos v = (X^2 - Y^2)/N, sin v = 2XY/N for the
rotated vector (X, Y) need no square root.

Run from the top-level FLINT directory:

    python3 dev/gen_mp_real_dense_tables.py
"""
import mpmath
from mpmath import mp, mpf

mp.prec = 1200

B = 2 ** 64
W = 10              # fraction limbs stored per value
BITS = 7            # bits per reduction level
EXP_LEVELS = 9      # G_2 .. G_EXP_LEVELS
LOG_LEVELS = 9      # M_1 .. M_LOG_LEVELS
ATAN_LEVELS = 6     # A_1 .. A_ATAN_LEVELS
EXP_TOP_J = 3       # the exp top level (after EXP_LEVELS): j <= EXP_TOP_J

LICENSE = """/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/"""


def limbs_of(x, integral):
    """floor(x B^W) as W fraction limbs (+ an integral limb)"""
    v = int(mpmath.floor(x * mpf(B) ** W))
    assert v >= 0
    out = []
    for _ in range(W):
        out.append(v % B)
        v //= B
    if integral:
        assert v < B
        out.append(v)
    else:
        assert v == 0, x
    return out


TABLES = []         # (name, comment, width, data)
DEFINES = []        # (name, value, comment)


def table(name, comment, values, integral=False):
    width = W + (1 if integral else 0)
    data = []
    for x in values:
        data.extend(limbs_of(x, integral))
    TABLES.append((name, comment, width, data))


def exp_index_bounds(levels):
    """J_k (the largest index at level k) for k = 2 .. levels, and the
    final remainder bound after each level"""
    # remainder before level 2: below 2^-7 + 2^-63 (L_1 one ulp low)
    rmax = mpf(2) ** -BITS + mpf(2) ** -63
    J = {}
    after = {}
    ulp = mpf(2) ** -(64 * 4)          # n >= 4 fraction limbs
    for k in range(2, levels + 1):
        s = mpf(2) ** (-BITS * k)
        jm = int(mpmath.floor(rmax / s))
        if jm * s == rmax:
            jm -= 1                     # rmax is exclusive
        J[k] = jm
        c = jm * s
        # t - G < (j + 1) s - log(1 + j s) + (k - 1) ulps, increasing in j
        rmax = s + (c - mpmath.log1p(c)) + k * ulp + mpf(2) ** -200
        after[k] = rmax
    return J, after


def build_exp():
    J, after = exp_index_bounds(EXP_LEVELS)
    c1 = []
    l1 = []
    s1 = min(j for j in range(128) if mpmath.exp(mpf(j) / 2 ** BITS) >= 2)
    for j in range(128):
        sh = 1 if j >= s1 else 0
        c = int(mpmath.floor((mpmath.exp(mpf(j) / 2 ** BITS) / 2 ** sh - 1) * B))
        assert 0 <= c < B
        c1.append(c)
        lv = sh * mpmath.log(2) + mpmath.log1p(mpf(c) / B)
        assert lv <= mpf(j) / 2 ** BITS < lv + mpf(2) ** -63
        l1.append(lv)
    DEFINES.append(("MP_REAL_DENSE_EXP_S1", s1,
                    "the first j with exp(j/2^7) >= 2"))
    TABLES.append(("_mp_real_dense_exp_c1",
        "C_1[j] = floor((exp(j/2^7) / 2^s - 1) 2^64), j < 128", 1, c1))
    # L_1 and G_2 are in dense_rt_tables.c (32 limbs, shared with the
    # run-time kernels, see dev/gen_mp_real_dense_rt_tables.py)
    values = []
    off = 0
    for k in range(2, EXP_LEVELS + 1):
        s = mpf(2) ** (-BITS * k)
        if k <= 2:
            assert J[k] <= 128
            continue
        DEFINES.append(("MP_REAL_DENSE_EXP_G%d" % k, off,
                        "limb offset of G_%d (j <= %d)" % (k, J[k])))
        for j in range(J[k] + 1):
            values.append(mpmath.log1p(j * s))
        off += (J[k] + 1) * W
        # the remainder after the last level must stay below 2^-(7k - 1)
        assert after[k] < mpf(2) ** -(BITS * k - 1)
    table("_mp_real_dense_exp_g",
          "log(1 + j 2^-7k), j <= J_k, k = 3 .. %d (offsets "
          "MP_REAL_DENSE_EXP_G<k>)" % EXP_LEVELS, values)
    DEFINES.append(("MP_REAL_DENSE_EXP_LEVELS", EXP_LEVELS,
                    "levels covered (the first one by C_1, L_1)"))
    # the top level after the last one: j = the whole top limb of t
    # (t < 2^-63 (1 + 2^-8), so j <= 2), t <- t - log(1 + j / 2^64),
    # leaving t below 2^-64 + 2^-127 + one ulp (below 2^-64 unless the
    # limb below the top is within a few units of 2^64)
    jt = int(mpmath.floor(after[EXP_LEVELS] * B))
    assert jt <= EXP_TOP_J
    table("_mp_real_dense_exp_gtop",
          "log(1 + j 2^-64), j <= %d (the top level after the last of "
          "seven bits)" % EXP_TOP_J,
          [mpmath.log1p(mpf(j) / B) for j in range(EXP_TOP_J + 1)])
    DEFINES.append(("MP_REAL_DENSE_EXP_TOP_J", EXP_TOP_J,
                    "the largest index of the top level"))


def log_d(k, j):
    """D_k(j) as an exact rational (one limb)"""
    if k == 1:
        return mpf((B * j) // (2 ** BITS + j)) / B
    s = BITS * k
    c2 = j * j * B                       # c^2 2^64 = j^2 2^64 / 2^(2s)
    num = -(-c2 // 2 ** (2 * s))         # ceil
    return mpf(j) / 2 ** s - mpf(num) / B


def build_log1p():
    ulp = mpf(2) ** -(64 * 4)
    values = []
    off = 0
    rmax = mpf(1)                        # u < 1 before level 1
    dvals = []
    for k in range(1, LOG_LEVELS + 1):
        s = mpf(2) ** (-BITS * k)
        jm = int(mpmath.floor(rmax / s))
        if jm * s == rmax:
            jm -= 1
        if k == 1:
            jm = 127
        # M_1 and M_2 are in dense_rt_tables.c (j <= 127 resp. 128)
        if k >= 3:
            DEFINES.append(("MP_REAL_DENSE_LOG_M%d" % k, off,
                            "limb offset of M_%d (j <= %d)" % (k, jm)))
        else:
            assert jm <= 128
        worst = mpf(0)
        for j in range(jm + 1):
            d = log_d(k, j)
            c = j * s
            assert 0 <= d <= c / (1 + c)
            if k >= 3:
                values.append(-mpmath.log1p(-d))
            if k == 1:
                dvals.append(d)
            worst = max(worst, (1 + (j + 1) * s) * (1 - d) - 1)
        if k >= 3:
            off += (jm + 1) * W
        rmax = worst + (k + 1) * ulp + mpf(2) ** -200
        # after the last level the series covers u < 2^-(7k - 1)
        if k >= 2:
            assert rmax < mpf(2) ** -(BITS * k - 1)
    TABLES.append(("_mp_real_dense_log_d1",
        "D_1[j] = floor(2^64 j / (2^7 + j)), j < 128", 1,
        [int(d * B) for d in dvals]))
    table("_mp_real_dense_log_m",
          "-log(1 - D_k(j)), j <= K_k, k = 3 .. %d (offsets "
          "MP_REAL_DENSE_LOG_M<k>)" % LOG_LEVELS, values)
    DEFINES.append(("MP_REAL_DENSE_LOG_LEVELS", LOG_LEVELS, "log1p levels"))


def build_atan():
    # A_1 and A_2 are in dense_rt_tables.c (j < 128 resp. j <= 133)
    values = []
    off = 0
    # entering level k >= 2, q < g 2^-7(k-1): the estimate of q is low
    # by at most rel (2^-12 deliberately, the drift of X after level 1,
    # truncations), so q_k < 2^-7k + rel q_{k-1}, with the fixed point
    # g = 1 / (1 - 2^7 rel); j_k <= floor(g 2^7)
    rel = mpf(2) ** -12 + mpf(2) ** -13.8 + mpf(2) ** -55
    g = 1 / (1 - 2 ** BITS * rel)
    jm = int(mpmath.floor(g * 2 ** BITS))
    assert g * 2 ** -(BITS * ATAN_LEVELS) < mpf(2) ** -(BITS * ATAN_LEVELS - 1)
    assert jm <= 133
    for k in range(3, ATAN_LEVELS + 1):
        DEFINES.append(("MP_REAL_DENSE_ATAN_A%d" % k, off,
                        "limb offset of A_%d (j <= %d)" % (k, jm)))
        for j in range(jm + 1):
            values.append(mpmath.atan(mpf(j) / 2 ** (BITS * k)))
        off += (jm + 1) * W
    table("_mp_real_dense_atan_a",
          "atan(j 2^-7k), k = 3 .. %d (offsets MP_REAL_DENSE_ATAN_A<k>)"
          % ATAN_LEVELS, values)
    DEFINES.append(("MP_REAL_DENSE_ATAN_LEVELS", ATAN_LEVELS, "atan levels"))


def build_sin_cos():
    c1 = []
    values = []
    for j in range(64):
        c = int(mpmath.floor(mpmath.tan(mpf(j) / 2 ** BITS) * B))
        assert c < B
        c1.append(c)
        a = mpmath.atan(mpf(c) / B)
        assert a <= mpf(j) / 2 ** BITS
        assert mpf(j) / 2 ** BITS - a < mpf(2) ** -63
        values.append(a)
    TABLES.append(("_mp_real_dense_sin_cos_c1",
        "C_1[j] = floor(tan(j/2^7) 2^64), j < 64", 1, c1))
    # T_1 is in dense_rt_tables.c


def emit():
    c = open("src/mp_real/dense_tables.c", "w")
    h = open("src/mp_real/dense_tables.h", "w")
    hdr = "/* GENERATED by dev/gen_mp_real_dense_tables.py -- do not edit */\n"
    total = sum(len(d) for _, _, _, d in TABLES) * 8
    c.write(LICENSE + "\n\n" + hdr + "\n")
    c.write('#include "flint.h"\n#include "dense_tables.h"\n\n'
            '#if FLINT_BITS == 64\n\n')
    h.write(LICENSE + "\n\n" + hdr + "\n")
    h.write("#ifndef MP_REAL_DENSE_TABLES_H\n#define MP_REAL_DENSE_TABLES_H\n\n")
    h.write('#include "flint.h"\n\n#if FLINT_BITS == 64\n\n')
    h.write("/* %d bytes in total */\n" % total)
    h.write("#define MP_REAL_DENSE_W %d\n" % W)
    for name, value, comment in DEFINES:
        h.write("/* %s */\n#define %s %d\n" % (comment, name, value))
    h.write("\n")
    for name, comment, width, data in TABLES:
        h.write("/* %s; %d limbs per entry */\n" % (comment, width))
        h.write("FLINT_DLL extern const ulong %s[%d];\n" % (name, len(data)))
        c.write("/* %s; %d limbs per entry */\n" % (comment, width))
        c.write("const ulong %s[%d] = {\n" % (name, len(data)))
        for i in range(0, len(data), width):
            c.write("    " + ", ".join("UWORD(0x%016x)" % d
                                       for d in data[i:i + width]) + ",\n")
        c.write("};\n\n")
    c.write("#endif\n")
    h.write("\n#endif\n\n#endif\n")


if __name__ == "__main__":
    build_exp()
    build_log1p()
    build_atan()
    build_sin_cos()
    emit()
