#!/usr/bin/env python3
"""Emit src/mp_real/elem_tables.c and src/mp_real/elem_tables.h: the
static tables of the elementary functions of the mp_real module, for
64-bit and 32-bit limbs.

  * floor(c B^N) for c = pi/4, log 2, 2/pi at N = CONST_N limbs (floors
    nest: the top n limbs are floor(c B^n) for n <= N);
  * the coefficients of the tapered series (series_tapered.c):

        tan:   c_j = T_(j+2) / (2j+3)!   (tangent numbers T_m)
        atan:  c_j = 1 / (2j+3)          (shared with atanh)
        sin:   c_j = 1 / (2j+3)!
        cos:   c_j = 1 / (2j+2)!          (1 - cos t = z V_0)

    for sizes n <= NMAX limbs and reduced arguments t < 2^-r, r >= RMIN:
    (the tangent's leading terms also on common denominators, in chunks
    of 1 .. TAN_RS_WMAX limbs, for series_rs.c; see emit_tan_rs)
    level j is read at L_j = ceil((64 n + g - r e_j) / 64) <= W_j limbs
    (e_j = 2j + 3, resp. 2j + 2 for cos; g = bits(K) + 1 at most), so
    coefficient j is stored as floor(c_j B^W_j) in W_j limbs, whose top
    L limbs are floor(c_j B^L); K is the number of levels at (NMAX,
    RMIN), the worst case; lg_j is an upper bound for log2 c_j.

Run from the top-level FLINT directory:

    python3 dev/gen_mp_real_elem_tables.py
"""
from fractions import Fraction
from math import factorial

CONST_N = {64: 64, 32: 128}

# (name, NMAX, RMIN) per family; NMAX in limbs of either size
FAMILIES = [("tan", 80, 32), ("atan", 13, 10), ("sin", 16, 18), ("cos", 16, 18)]

# chunk widths 1 .. TAN_RS_WMAX limbs for the tangent's common-denominator
# rectangular splitting
TAN_RS_WMAX = 12

LICENSE = """/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/"""


def tan_coeffs(K):
    # c_j = coefficient of t^(2j+3) in tan t, via tan' = 1 + tan^2
    a = [Fraction(1)]                      # a_m: coefficient of t^(2m+1)
    for m in range(1, K + 2):
        a.append(sum(a[i] * a[m - 1 - i] for i in range(m)) / (2 * m + 1))
    return [a[j + 1] for j in range(K)]


def coeffs(fam, K):
    if fam == "tan":
        return tan_coeffs(K)
    if fam == "atan":
        return [Fraction(1, 2 * j + 3) for j in range(K)]
    if fam == "sin":
        return [Fraction(1, factorial(2 * j + 3)) for j in range(K)]
    return [Fraction(1, factorial(2 * j + 2)) for j in range(K)]


def lg_upper(c):
    """smallest integer l with c < 2^l"""
    l = c.numerator.bit_length() - c.denominator.bit_length()
    while Fraction(2) ** l <= c:
        l += 1
    while Fraction(2) ** (l - 1) > c:
        l -= 1
    return l


def emit_family(o, h, fam, nmax, rmin, bits):
    e0 = 2 if fam == "cos" else 3
    # levels at the worst case (nmax, rmin), as series_tapered.c counts them
    cs = coeffs(fam, 400)
    K = 0
    while lg_upper(cs[K]) - rmin * (e0 + 2 * K) >= -bits * nmax - 2:
        K += 1
    K = max(K, 1)
    g = K.bit_length() + 1
    W = []
    for j in range(K):
        b = bits * nmax + g - rmin * (e0 + 2 * j)
        W.append(1 if b <= 0 else min(nmax, -(-b // bits)))
    for j in range(K - 1, 0, -1):
        W[j - 1] = max(W[j - 1], W[j])
    off = []
    flat = []
    for j in range(K):
        off.append(len(flat))
        v = (cs[j].numerator << (bits * W[j])) // cs[j].denominator
        flat += [(v >> (bits * i)) & ((1 << bits) - 1) for i in range(W[j])]
    fmt = "UWORD(0x%016x)" if bits == 64 else "UWORD(0x%08x)"
    up = fam.upper()
    h("#define MP_REAL_SERIES_%s_NMAX %d" % (up, nmax))
    h("#define MP_REAL_SERIES_%s_RMIN %d" % (up, rmin))
    h("#define MP_REAL_SERIES_%s_K %d" % (up, K))
    h("FLINT_DLL extern const ulong _mp_real_series_%s_c[%d];" % (fam, len(flat)))
    h("FLINT_DLL extern const short _mp_real_series_%s_off[%d];" % (fam, K))
    h("FLINT_DLL extern const short _mp_real_series_%s_w[%d];" % (fam, K))
    h("FLINT_DLL extern const short _mp_real_series_%s_lg[%d];" % (fam, K))
    o("const ulong _mp_real_series_%s_c[%d] = {" % (fam, len(flat)))
    for i in range(0, len(flat), 4 if bits == 64 else 6):
        o("    " + ", ".join(fmt % x for x in flat[i:i + (4 if bits == 64 else 6)]) + ",")
    o("};")
    o("const short _mp_real_series_%s_off[%d] = {%s};" % (fam, K, ", ".join(map(str, off))))
    o("const short _mp_real_series_%s_w[%d] = {%s};" % (fam, K, ", ".join(map(str, W))))
    o("const short _mp_real_series_%s_lg[%d] = {%s};" % (fam, K, ", ".join(str(lg_upper(c)) for c in cs[:K])))
    o("")
    return len(flat)


def emit_tan_rs(o, h, bits):
    """chunks of the tangent series on common denominators: chunk j
    (w = j + 1 limbs) holds the terms k in [start_j, start_(j+1)) whose
    cumulative lcm Q_j of the reduced denominators of c_0 .. c_k stays
    below 2^(bits w - 2); numerators N_k = c_k Q_j (w limbs), Q_j and
    R_j = Q_j / Q_(j-1) (w limbs, R_0 = Q_0)"""
    from math import lcm
    cs = tan_coeffs(200)
    starts = [0]
    Qs = []
    Q = 1
    k = 0
    for j in range(TAN_RS_WMAX):
        w = j + 1
        while k < len(cs) and lcm(Q, cs[k].denominator) < (1 << (bits * w - 2)):
            Q = lcm(Q, cs[k].denominator)
            k += 1
        Qs.append(Q)
        starts.append(k)
    J = len(Qs)
    mask = (1 << bits) - 1
    fmt = "UWORD(0x%016x)" if bits == 64 else "UWORD(0x%08x)"
    def limbs(v, w):
        return [(v >> (bits * i)) & mask for i in range(w)]
    qflat, rflat, nflat, noff = [], [], [], []
    for j in range(J):
        w = j + 1
        qflat += limbs(Qs[j], w)
        R = Qs[j] // (Qs[j - 1] if j else 1)
        assert R * (Qs[j - 1] if j else 1) == Qs[j]
        rflat += limbs(R, w)
        for kk in range(starts[j], starts[j + 1]):
            v = cs[kk] * Qs[j]
            assert v.denominator == 1 and v.numerator < (1 << (bits * w))
            noff.append(len(nflat))
            nflat += limbs(v.numerator, w)
    h("#define MP_REAL_SERIES_TAN_RS_J %d" % J)
    h("FLINT_DLL extern const short _mp_real_series_tan_rs_start[%d];" % (J + 1))
    h("FLINT_DLL extern const ulong _mp_real_series_tan_rs_q[%d];" % len(qflat))
    h("FLINT_DLL extern const ulong _mp_real_series_tan_rs_r[%d];" % len(rflat))
    h("FLINT_DLL extern const ulong _mp_real_series_tan_rs_num[%d];" % len(nflat))
    h("FLINT_DLL extern const short _mp_real_series_tan_rs_num_off[%d];" % len(noff))
    o("/* tangent chunks on common denominators (series_rs.c) */")
    o("const short _mp_real_series_tan_rs_start[%d] = {%s};" % (J + 1, ", ".join(map(str, starts))))
    for name, arr in (("q", qflat), ("r", rflat), ("num", nflat)):
        o("const ulong _mp_real_series_tan_rs_%s[%d] = {" % (name, len(arr)))
        step = 4 if bits == 64 else 6
        for i in range(0, len(arr), step):
            o("    " + ", ".join(fmt % x for x in arr[i:i + step]) + ",")
        o("};")
    o("const short _mp_real_series_tan_rs_num_off[%d] = {%s};" % (len(noff), ", ".join(map(str, noff))))
    o("")


def const_limbs(c, bits, N):
    v = int(c * (1 << (bits * N)))
    return [(v >> (bits * i)) & ((1 << bits) - 1) for i in range(N)]


def main():
    from mpmath import mp, pi, log
    out = []
    hdr = []
    o = out.append
    h = hdr.append
    o(LICENSE)
    o("")
    o("/* GENERATED by dev/gen_mp_real_elem_tables.py -- do not edit */")
    o("")
    o('#include "flint.h"')
    o('#include "mp_real.h"')
    o('#include "impl.h"')
    o('#include "elem_tables.h"')
    o("")
    h(LICENSE)
    h("")
    h("/* GENERATED by dev/gen_mp_real_elem_tables.py -- do not edit */")
    h("")
    h("#ifndef MP_REAL_ELEM_TABLES_H")
    h("#define MP_REAL_ELEM_TABLES_H")
    h("")
    for bits in (64, 32):
        o("#if FLINT_BITS == %d" % bits)
        o("")
        h("#if FLINT_BITS == %d" % bits)
        N = CONST_N[bits]
        h("#define MP_REAL_CONST_STATIC_N %d" % N)
        mp.prec = bits * N + 256
        fmt = "UWORD(0x%016x)" if bits == 64 else "UWORD(0x%08x)"
        for name, val in (("pi4", pi / 4), ("log2", log(2)), ("2_div_pi", 2 / pi)):
            L = const_limbs(val, bits, N)
            h("FLINT_DLL extern const ulong _mp_real_const_%s_static[%d];" % (name, N))
            o("/* floor(%s B^%d), low limb first */" % (name, N))
            o("const ulong _mp_real_const_%s_static[%d] = {" % (name, N))
            step = 4 if bits == 64 else 6
            for i in range(0, N, step):
                o("    " + ", ".join(fmt % x for x in L[i:i + step]) + ",")
            o("};")
            o("")
        for fam, nmax, rmin in FAMILIES:
            emit_family(o, h, fam, nmax, rmin, bits)
        emit_tan_rs(o, h, bits)
        o("#endif")
        o("")
        h("#endif")
        h("")
    h("#endif")
    open("src/mp_real/elem_tables.c", "w").write("\n".join(out) + "\n")
    open("src/mp_real/elem_tables.h", "w").write("\n".join(hdr) + "\n")


if __name__ == "__main__":
    main()
