/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "arf.h"

/*
    Complex reciprocals, quotients and square roots on mantissa/exponent
    pairs, in the _approx flavour (rn limbs per part, each within about
    an ulp of the result's magnitude, the parts truncated): the
    compositions of nfloat_complex on the kernels, every intermediate
    quantity held at sn = rn + 1 limbs (the guard limb keeps the few ulps
    of the compositions below the ulp of the output; without it the parts
    were within a few ulps, which acf's callers, at 2^(2 - prec), do not
    allow; the guard limb costs 30-150% at small sizes against the
    compositions at rn limbs, but the kernels remain well ahead of the
    general arf arithmetic they replace). The exponents are slongs, far
    from overflow (the callers scale extreme operands first).

      1 / (a + b i) = (a - b i) t,  t = 1 / s,  s = a^2 + b^2 (one
                      reciprocal and two products: a division costs
                      several products);
      (a + b i) / (c + d i) = (u + v i) t,  u = a c + b d,
                              v = b c - a d,  s = c^2 + d^2;
      sqrt(a + b i) = (q, b / (2 q)) for a >= 0,  q = sqrt((r + |a|) / 2),
                      r = |a + b i|, and (|b| / (2 q), sgn(b) q) for
                      a < 0 (the larger part first: no cancellation);
      1 / sqrt(a + b i) = (T v, -b v) for a >= 0, T = r + |a|,
                          v = 1 / sqrt(2 r^2 T) = 1 / sqrt(r (T^2 + b^2)),
                          and (|b| v, -sgn(b) T v) for a < 0.

    The functions write the magnitudes of the parts; the signs follow
    the rules above and are the caller's (the imaginary part of the
    square root takes the sign of b, that of the reciprocal root the
    opposite).
*/

/* The sums of the compositions, truncated at rn limbs: for magnitudes
   X = {x, xn} B^-xn and Y = {y, yn} B^-yn in [1/2, 1) with Y scaled by
   2^-shift (shift >= 0), z = X + Y 2^-shift (_add_trunc) or
   X - Y 2^-shift (_sub_trunc, with X >= Y 2^-shift unless shift = 0,
   in which case the sign of the difference comes out through *neg),
   within one unit of the rn-limb window of the result (the bits of Y
   below the window are dropped, and borrow from them when X has none
   there). The window of Y lands in the scratch t of rn limbs. *zexp is
   the correction of the exponent of X (WORD_MIN for an exact zero). */

/* the window [0, rn) of Y 2^-shift into t (which may be y itself): the
   top rn - shift / B limbs of Y, shifted right by shift mod B, at the
   bottom of t; returns whether any bit of Y below the window is set */
static int
_shifted_window(nn_ptr t, nn_srcptr y, slong yn, slong shift, slong rn)
{
    slong ls = shift / FLINT_BITS, m = rn - ls, i;
    unsigned int b = shift % FLINT_BITS;
    int sticky = 0;

    if (m <= 0)
    {
        _flint_mpn_zero_small(t, rn);
        return 1;
    }

    /* the top m limbs of Y (zero padded below) into t[0, m) */
    if (yn >= m)
    {
        for (i = 0; i < yn - m; i++)
            sticky |= (y[i] != 0);
        memmove(t, y + yn - m, m * sizeof(ulong));
    }
    else
    {
        memmove(t + m - yn, y, yn * sizeof(ulong));
        _flint_mpn_zero_small(t, m - yn);
    }
    _flint_mpn_zero_small(t + m, ls);

    if (b != 0)
        sticky |= (mpn_rshift(t, t, m, b) != 0);

    return sticky;
}

/* the top rn limbs of X into z (zero padded below when xn < rn) */
static void
_top_window(nn_ptr z, nn_srcptr x, slong xn, slong rn)
{
    if (xn >= rn)
        flint_mpn_copyi(z, x + xn - rn, rn);
    else
    {
        flint_mpn_copyi(z + rn - xn, x, xn);
        _flint_mpn_zero_small(z, rn - xn);
    }
}

static void
_add_trunc(nn_ptr z, slong * zexp, nn_ptr t, nn_srcptr x, slong xn, nn_srcptr y, slong yn, slong shift, slong rn)
{
    _top_window(z, x, xn, rn);
    _shifted_window(t, y, yn, shift, rn);

    if (mpn_add_n(z, z, t, rn))
    {
        mpn_rshift(z, z, rn, 1);
        z[rn - 1] |= LIMB_TOP;
        *zexp = 1;
    }
    else
        *zexp = 0;
}

static void
_sub_trunc(nn_ptr z, slong * zexp, int * neg, nn_ptr t, nn_srcptr x, slong xn, nn_srcptr y, slong yn, slong shift, slong rn)
{
    unsigned int lead;
    int sticky;

    _top_window(z, x, xn, rn);
    sticky = _shifted_window(t, y, yn, shift, rn) && (xn <= rn);
    *neg = 0;

    if (shift == 0 && mpn_cmp(z, t, rn) < 0)
    {
        mpn_sub_n(z, t, z, rn);
        *neg = 1;
        sticky = 0;
    }
    else
        mpn_sub_n(z, z, t, rn);

    if (sticky)
        mpn_sub_1(z, z, rn, 1);

    if (z[rn - 1] == 0)
    {
        slong k = rn - 1;
        while (k >= 0 && z[k] == 0)
            k--;
        if (k < 0)
        {
            *zexp = WORD_MIN;
            return;
        }
        /* k + 1 nonzero limbs: moved up */
        flint_mpn_copyd(z + rn - k - 1, z, k + 1);
        _flint_mpn_zero_small(z, rn - k - 1);
        lead = flint_clz(z[rn - 1]) + (rn - k - 1) * FLINT_BITS;
        if (lead % FLINT_BITS != 0)
            mpn_lshift(z, z, rn, lead % FLINT_BITS);
        *zexp = -(slong) lead;
        return;
    }

    lead = flint_clz(z[rn - 1]);
    if (lead != 0)
        mpn_lshift(z, z, rn, lead);
    *zexp = -(slong) lead;
}

/* z = (-1)^psgn p + (-1)^qsgn q for magnitudes p, q (rn limbs each) with
   exponents, at rn limbs (truncated): the sum or the difference of the
   magnitudes, the larger exponent first; *zneg the sign, *zexp = WORD_MIN
   for an exact zero; p and q lie in scratch (the smaller is clobbered) */
static void
_sum2(nn_ptr z, slong * zexp, int * zneg, nn_ptr p, slong pexp, int psgn, nn_ptr q, slong qexp, int qsgn, slong rn)
{
    slong fix;
    int neg;

    if (pexp < qexp)
    {
        nn_ptr tp = p; slong te = pexp; int ts = psgn;
        p = q; pexp = qexp; psgn = qsgn;
        q = tp; qexp = te; qsgn = ts;
    }

    if (psgn == qsgn)
    {
        _add_trunc(z, &fix, q, p, rn, q, rn, pexp - qexp, rn);
        *zneg = psgn;
    }
    else
    {
        _sub_trunc(z, &fix, &neg, q, p, rn, q, rn, pexp - qexp, rn);
        *zneg = psgn ^ neg;
    }
    *zexp = (fix == WORD_MIN) ? WORD_MIN : pexp + fix;
}

/* s = a^2 + b^2 at rn limbs (truncated) with its exponent; scratch t of
   2 rn limbs */
FLINT_FORCE_INLINE void
_sosq(nn_ptr s, slong * sexp, nn_ptr t, nn_srcptr a, slong an, slong ea, nn_srcptr b, slong bn, slong eb, slong rn)
{
    slong pexp, qexp;
    int neg;
    _arf_mpn_mul_approx(t, &pexp, a, an, a, an, rn * FLINT_BITS);
    _arf_mpn_mul_approx(t + rn, &qexp, b, bn, b, bn, rn * FLINT_BITS);
    _sum2(s, sexp, &neg, t, 2 * ea + pexp, 0, t + rn, 2 * eb + qexp, 0, rn);
}

/* z = x y truncated at rn limbs through scratch t (z may alias an
   operand, which the product kernels do not allow) */
FLINT_FORCE_INLINE void
_mul_via(nn_ptr z, slong * zexp, nn_ptr t, nn_srcptr x, slong xn, slong ex, nn_srcptr y, slong yn, slong ey, slong rn)
{
    slong fix;
    _arf_mpn_mul_approx(t, &fix, x, xn, y, yn, rn * FLINT_BITS);
    flint_mpn_copyi(z, t, rn);
    *zexp = ex + ey + fix;
}

void
_arf_mpn_cinv_approx(nn_ptr re, slong * reexp, nn_ptr im, slong * imexp,
    nn_srcptr a, slong an, slong ea, nn_srcptr b, slong bn, slong eb, slong prec)
{
    slong rn = (prec + FLINT_BITS - 1) / FLINT_BITS, sn = rn + 1, sexp, texp, fix;
    nn_ptr s, t, w;
    TMP_INIT;

    /* (the scratch of the kernels called below is not ours to use) */
    TMP_START;
    s = TMP_ALLOC(4 * sn * sizeof(ulong));
    t = s + sn;
    w = t + sn;

    _sosq(s, &sexp, w, a, an, ea, b, bn, eb, sn);
    _arf_mpn_inv_approx(t, &fix, s, sn, sn * FLINT_BITS);
    texp = fix - sexp;

    _mul_via(re, reexp, s, a, an, ea, t, sn, texp, rn);
    _mul_via(im, imexp, s, b, bn, eb, t, sn, texp, rn);
    TMP_END;
}

void
_arf_mpn_cdiv_approx(nn_ptr re, slong * reexp, int * reneg, nn_ptr im, slong * imexp, int * imneg,
    nn_srcptr a, slong an, slong ea, int sa, nn_srcptr b, slong bn, slong eb, int sb,
    nn_srcptr c, slong cn, slong ec, int sc, nn_srcptr d, slong dn, slong ed, int sd, slong prec)
{
    slong rn = (prec + FLINT_BITS - 1) / FLINT_BITS, sn = rn + 1, sexp, texp, uexp, vexp, pexp, qexp, fix;
    nn_ptr s, t, u, v, w;
    int uneg, vneg;
    TMP_INIT;

    TMP_START;
    s = TMP_ALLOC(6 * sn * sizeof(ulong));
    t = s + sn;
    u = t + sn;
    v = u + sn;
    w = v + sn;

    /* u = a c + b d, v = b c - a d, with their signs */
    _arf_mpn_mul_approx(w, &pexp, a, an, c, cn, sn * FLINT_BITS);
    _arf_mpn_mul_approx(w + sn, &qexp, b, bn, d, dn, sn * FLINT_BITS);
    _sum2(u, &uexp, &uneg, w, ea + ec + pexp, sa ^ sc, w + sn, eb + ed + qexp, sb ^ sd, sn);
    _arf_mpn_mul_approx(w, &pexp, b, bn, c, cn, sn * FLINT_BITS);
    _arf_mpn_mul_approx(w + sn, &qexp, a, an, d, dn, sn * FLINT_BITS);
    _sum2(v, &vexp, &vneg, w, eb + ec + pexp, sb ^ sc, w + sn, ea + ed + qexp, !(sa ^ sd), sn);

    _sosq(s, &sexp, w, c, cn, ec, d, dn, ed, sn);
    _arf_mpn_inv_approx(t, &fix, s, sn, sn * FLINT_BITS);
    texp = fix - sexp;

    if (uexp == WORD_MIN)
        *reexp = WORD_MIN;
    else
        _mul_via(re, reexp, s, u, sn, uexp, t, sn, texp, rn);
    *reneg = uneg;

    if (vexp == WORD_MIN)
        *imexp = WORD_MIN;
    else
        _mul_via(im, imexp, s, v, sn, vexp, t, sn, texp, rn);
    *imneg = vneg;
    TMP_END;
}

/* r = |a + b i| at rn limbs: the root of a^2 + b^2; scratch s of 3 rn
   limbs */
FLINT_FORCE_INLINE void
_cabs(nn_ptr r, slong * rexp, nn_ptr s, nn_srcptr a, slong an, slong ea, nn_srcptr b, slong bn, slong eb, slong rn)
{
    slong sexp, fix;
    int t;
    _sosq(s, &sexp, s + rn, a, an, ea, b, bn, eb, rn);
    t = sexp & 1;
    _arf_mpn_sqrt_approx(r, &fix, s, rn, t, rn * FLINT_BITS);
    *rexp = ((sexp + t) >> 1) + fix;
}

/* T = r + |a| at sn limbs (r >= |a|, but the truncated r may fall just
   below a power of two |a|: the operand with the larger exponent first);
   scratch t of sn limbs */
FLINT_FORCE_INLINE void
_cabs_plus(nn_ptr T, slong * Texp, nn_ptr t, nn_srcptr r, slong rexp, slong sn, nn_srcptr a, slong an, slong ea)
{
    slong fix;
    if (rexp >= ea)
    {
        _add_trunc(T, &fix, t, r, sn, a, an, rexp - ea, sn);
        *Texp = rexp + fix;
    }
    else
    {
        _add_trunc(T, &fix, t, a, an, r, sn, ea - rexp, sn);
        *Texp = ea + fix;
    }
}

void
_arf_mpn_csqrt_approx(nn_ptr re, slong * reexp, nn_ptr im, slong * imexp,
    nn_srcptr a, slong an, slong ea, int sa, nn_srcptr b, slong bn, slong eb, slong prec)
{
    slong rn = (prec + FLINT_BITS - 1) / FLINT_BITS, sn = rn + 1;
    slong rexp, Texp, qexp, fix;
    nn_ptr s, r, T, q, big, small;
    slong * bigexp, * smallexp;
    int t;
    TMP_INIT;

    TMP_START;
    s = TMP_ALLOC(6 * sn * sizeof(ulong));
    r = s + 3 * sn;
    T = r + sn;
    q = T + sn;

    _cabs(r, &rexp, s, a, an, ea, b, bn, eb, sn);
    _cabs_plus(T, &Texp, s, r, rexp, sn, a, an, ea);

    /* q = sqrt(T / 2), at sn limbs for the quotient below, its top rn
       limbs the larger part */
    Texp -= 1;
    t = Texp & 1;
    _arf_mpn_sqrt_approx(q, &fix, T, sn, t, sn * FLINT_BITS);
    qexp = ((Texp + t) >> 1) + fix;

    /* the larger part is q, the other |b| / (2 q) */
    if (!sa)
    {
        big = re; bigexp = reexp;
        small = im; smallexp = imexp;
    }
    else
    {
        big = im; bigexp = imexp;
        small = re; smallexp = reexp;
    }

    _arf_mpn_div_approx(small, &fix, b, bn, q, sn, prec);
    *smallexp = eb - qexp - 1 + fix;
    flint_mpn_copyi(big, q + 1, rn);
    *bigexp = qexp;
    TMP_END;
}

void
_arf_mpn_crsqrt_approx(nn_ptr re, slong * reexp, nn_ptr im, slong * imexp,
    nn_srcptr a, slong an, slong ea, int sa, nn_srcptr b, slong bn, slong eb, slong prec)
{
    slong rn = (prec + FLINT_BITS - 1) / FLINT_BITS, sn = rn + 1;
    slong rexp, Texp, wexp, vexp, fix;
    nn_ptr s, r, T, w, v, big, small;
    slong * bigexp, * smallexp;
    int t;
    TMP_INIT;

    TMP_START;
    s = TMP_ALLOC(7 * sn * sizeof(ulong));
    r = s + 3 * sn;
    T = r + sn;
    w = T + sn;
    v = w + sn;

    _cabs(r, &rexp, s, a, an, ea, b, bn, eb, sn);
    _cabs_plus(T, &Texp, s, r, rexp, sn, a, an, ea);

    /* w = T^2 + b^2, then w = r w = 2 r^2 T */
    _sosq(w, &wexp, s, T, sn, Texp, b, bn, eb, sn);
    _arf_mpn_mul_approx(s, &fix, w, sn, r, sn, sn * FLINT_BITS);
    wexp += rexp + fix;

    /* v = 1 / sqrt(w) */
    t = wexp & 1;
    _arf_mpn_rsqrt_approx(v, &fix, s, sn, t, sn * FLINT_BITS);
    vexp = (1 - t) - ((wexp - t) >> 1) + fix;

    if (!sa)
    {
        big = re; bigexp = reexp;
        small = im; smallexp = imexp;
    }
    else
    {
        big = im; bigexp = imexp;
        small = re; smallexp = reexp;
    }

    _mul_via(small, smallexp, s, b, bn, eb, v, sn, vexp, rn);
    _mul_via(big, bigexp, s, T, sn, Texp, v, sn, vexp, rn);
    TMP_END;
}
