/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef MPN_EXTRAS_FP_SMALL_H
#define MPN_EXTRAS_FP_SMALL_H

/* Register kernels for floating-point quotients, square roots and
   reciprocal square roots at one to three limbs, shared by the
   floating-point types (divhigh.c, sqrthigh.c, arf, nfloat). The
   root kernels are for 64-bit limbs only. The fractions are A = a / B^k
   in [1/2, 1) (top bit set) with t = 0 or 1 the parity of the exponent,
   so that the quantities sqrt(A / 2^t) and M = (2^(2 - t) A)^(-1/2) lie
   in [1/2, 1) resp. (1/2, 1]. */

#include <math.h>
#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"

/*
    One-limb floating-point division: q = floor(a 2^(FLINT_BITS - s) / b)
    for normalised limbs a and b, s = [a >= b], exactly (the top bit of q
    set); *rem receives the remainder (nonzero iff inexact). Hardware
    division where it is fast (FLINT_PREINVERT_LIMB_USE_NATIVE), otherwise
    a precomputed inverse (binv = FLINT_MPN_INVERT_LIMB).
*/
FLINT_FORCE_INLINE ulong
_flint_mpn_fdiv_1_kernel_preinv(ulong * rem, ulong * s, ulong a, ulong b, ulong binv)
{
    ulong t, nh, nl, q, r;

    t = (a >= b);
    nh = a >> t;
    nl = (a << (FLINT_BITS - 1)) & (-t);
    FLINT_MPN_UDIV_QR_2BY1(q, r, nh, nl, b, binv);
    *rem = r;
    *s = t;
    return q;
}

FLINT_FORCE_INLINE ulong
_flint_mpn_fdiv_1_kernel(ulong * rem, ulong * s, ulong a, ulong b)
{
#if FLINT_PREINVERT_LIMB_USE_NATIVE
    ulong t, nh, nl, q, r;

    t = (a >= b);
    nh = a >> t;
    nl = (a << (FLINT_BITS - 1)) & (-t);
    udiv_qrnnd(q, r, nh, nl, b);
    *rem = r;
    *s = t;
    return q;
#else
    ulong binv;
    FLINT_MPN_INVERT_LIMB(binv, b);
    return _flint_mpn_fdiv_1_kernel_preinv(rem, s, a, b, binv);
#endif
}

/*
    Two-limb floating-point division: with the fractions (a1, a0) / B^2
    and (b1, b0) / B^2 in [1/2, 1) and s = [a >= b], (q1, q0) =
    floor(a B^2 / (2^s b)) exactly, which has the top bit set. The
    dividend is shifted right by s bits (branch-free) before the two
    3/2 divisions; the remainder gives the inexactness.
*/
#define _FLINT_MPN_FDIV_2_NUMERATOR(n2, n1, n0, s, a1, a0, b1, b0) \
    do { \
        (s) = ((a1) > (b1)) | (((a1) == (b1)) & ((a0) >= (b0))); \
        (n2) = (a1) >> (s); \
        (n1) = ((a0) >> (s)) | (((a1) << (FLINT_BITS - 1)) & (-(s))); \
        (n0) = ((a0) << (FLINT_BITS - 1)) & (-(s)); \
    } while (0)

/* (q1, q0) = floor((a1, a0, 0, 0) / (2^s (b1, b0))); returns nonzero iff
   inexact. */
FLINT_FORCE_INLINE int
_flint_mpn_fdiv_2_kernel_preinv(ulong * q1p, ulong * q0p, ulong * sp, ulong a1, ulong a0, ulong b1, ulong b0, ulong dinv)
{
    ulong s, n2, n1, n0, q1, q0, r1, r0, t1, t0;

    _FLINT_MPN_FDIV_2_NUMERATOR(n2, n1, n0, s, a1, a0, b1, b0);
    FLINT_MPN_UDIV_QR_3BY2(q1, r1, r0, n2, n1, n0, b1, b0, dinv);
    FLINT_MPN_UDIV_QR_3BY2(q0, t1, t0, r1, r0, UWORD(0), b1, b0, dinv);
    *q1p = q1;
    *q0p = q0;
    *sp = s;
    return (t1 | t0) != 0;
}

FLINT_FORCE_INLINE int
_flint_mpn_fdiv_2_kernel(ulong * q1p, ulong * q0p, ulong * sp, ulong a1, ulong a0, ulong b1, ulong b0)
{
#if FLINT_PREINVERT_LIMB_USE_NATIVE
    ulong s, n2, n1, n0, q1, q0, r1, r0, t1, t0;

    _FLINT_MPN_FDIV_2_NUMERATOR(n2, n1, n0, s, a1, a0, b1, b0);
    FLINT_MPN_UDIV_QR_3BY2_HW(q1, r1, r0, n2, n1, n0, b1, b0);
    FLINT_MPN_UDIV_QR_3BY2_HW(q0, t1, t0, r1, r0, UWORD(0), b1, b0);
    *q1p = q1;
    *q0p = q0;
    *sp = s;
    return (t1 | t0) != 0;
#else
    return _flint_mpn_fdiv_2_kernel_preinv(q1p, q0p, sp, a1, a0, b1, b0, flint_mpn_preinv1(b1, b0));
#endif
}

/*
    Division by a three-limb divisor: with the fractions A = (a1, a0) / B^2
    and D = (d2, d1, d0) / B^3 in [1/2, 1) and s = [A >= D],
    (q2, q1, q0) = floor(A B^3 / (2^s D)) exactly, which has the top bit
    set; returns nonzero iff inexact. The dividend (a1, a0, 0, 0, 0, 0) is
    shifted right by s bits and divided by the schoolbook method, one
    3/2 division per quotient limb with the full remainder carried along
    (as _flint_mpn_tdiv_qr_small, with the sizes fixed). In the rare
    case of a partial remainder sharing its top two limbs with the
    divisor, the quotient limb is B - 1.
*/
#define _FLINT_MPN_FDIV_3_STEP(q, r2, r1, r0, x, d2, d1, d0, DIV3BY2) \
    do { \
        ulong __fd3_q, __fd3_u1, __fd3_u0, __fd3_ph, __fd3_pl, __fd3_t2, __fd3_t1, __fd3_t0, __fd3_c; \
        if (FLINT_UNLIKELY((r2) == (d2) && (r1) == (d1))) \
        { \
            /* q = B - 1: the new remainder (R B + x) - (B - 1) D = \
               D + x - (d0 - r0) B, which is in [0, D) */ \
            __fd3_q = ~UWORD(0); \
            __fd3_t0 = (d0) + (x); \
            __fd3_c = __fd3_t0 < (x); \
            __fd3_t1 = (d1) + __fd3_c; \
            __fd3_c = __fd3_t1 < __fd3_c; \
            __fd3_t2 = (d2) + __fd3_c; \
            __fd3_c = __fd3_t2 < __fd3_c; \
            __fd3_u1 = (d0) - (r0); \
            __fd3_u0 = __fd3_t1 < __fd3_u1; \
            __fd3_t1 -= __fd3_u1; \
            FLINT_ASSERT(__fd3_c == (__fd3_t2 < __fd3_u0)); \
            __fd3_t2 -= __fd3_u0; \
            (r2) = __fd3_t2; (r1) = __fd3_t1; (r0) = __fd3_t0; \
        } \
        else \
        { \
            DIV3BY2(__fd3_q, __fd3_u1, __fd3_u0, (r2), (r1), (r0), (d2), (d1)); \
            /* (u1, u0, x) - q d0 */ \
            umul_ppmm(__fd3_ph, __fd3_pl, __fd3_q, (d0)); \
            __fd3_t0 = (x) - __fd3_pl; \
            __fd3_c = (x) < __fd3_pl; \
            __fd3_t1 = __fd3_u0 - __fd3_ph; \
            __fd3_c = (__fd3_u0 < __fd3_ph) | (__fd3_t1 < __fd3_c); \
            __fd3_t1 -= ((x) < __fd3_pl); \
            __fd3_t2 = __fd3_u1 - __fd3_c; \
            if (__fd3_u1 < __fd3_c) \
            { \
                /* negative: add D back once (the 3/2 estimate is at \
                   most one too large) */ \
                __fd3_q--; \
                __fd3_t0 += (d0); \
                __fd3_c = __fd3_t0 < (d0); \
                __fd3_t1 += (d1) + __fd3_c; \
                __fd3_c = (__fd3_t1 < (d1)) | (__fd3_c & (__fd3_t1 == (d1))); \
                __fd3_t2 += (d2) + __fd3_c; \
            } \
            (r2) = __fd3_t2; (r1) = __fd3_t1; (r0) = __fd3_t0; \
        } \
        (q) = __fd3_q; \
    } while (0)

/* s = [(a1, a0, 0) >= (d2, d1, d0)] and the top three limbs of the
   dividend shifted right by s (below D); the limbs below are zero */
#define _FLINT_MPN_FDIV_3_NUMERATOR(r2, r1, r0, s, a1, a0, d2, d1, d0) \
    do { \
        (s) = ((a1) > (d2)) | (((a1) == (d2)) & (((a0) > (d1)) | (((a0) == (d1)) & ((d0) == 0)))); \
        (r2) = (a1) >> (s); \
        (r1) = ((a0) >> (s)) | (((a1) << (FLINT_BITS - 1)) & (-(s))); \
        (r0) = ((a0) << (FLINT_BITS - 1)) & (-(s)); \
    } while (0)

/* with dinv = flint_mpn_preinv1(d2, d1), for several divisions by D */
FLINT_FORCE_INLINE int
_flint_mpn_fdiv_3_kernel_preinv(ulong * q2p, ulong * q1p, ulong * q0p, ulong * sp,
    ulong a1, ulong a0, ulong d2, ulong d1, ulong d0, ulong dinv)
{
    ulong s, r2, r1, r0, x = 0, q2, q1, q0;
#define _FDIV_3_DIV3BY2(q, r1, r0, n2, n1, n0, d1, d0) FLINT_MPN_UDIV_QR_3BY2(q, r1, r0, n2, n1, n0, d1, d0, dinv)
    _FLINT_MPN_FDIV_3_NUMERATOR(r2, r1, r0, s, a1, a0, d2, d1, d0);
    _FLINT_MPN_FDIV_3_STEP(q2, r2, r1, r0, x, d2, d1, d0, _FDIV_3_DIV3BY2);
    _FLINT_MPN_FDIV_3_STEP(q1, r2, r1, r0, x, d2, d1, d0, _FDIV_3_DIV3BY2);
    _FLINT_MPN_FDIV_3_STEP(q0, r2, r1, r0, x, d2, d1, d0, _FDIV_3_DIV3BY2);
#undef _FDIV_3_DIV3BY2
    *q2p = q2;
    *q1p = q1;
    *q0p = q0;
    *sp = s;
    return (r2 | r1 | r0) != 0;
}

FLINT_FORCE_INLINE int
_flint_mpn_fdiv_3_kernel(ulong * q2p, ulong * q1p, ulong * q0p, ulong * sp,
    ulong a1, ulong a0, ulong d2, ulong d1, ulong d0)
{
#if FLINT_PREINVERT_LIMB_USE_NATIVE
    ulong s, r2, r1, r0, x = 0, q2, q1, q0;
#define _FDIV_3_DIV3BY2(q, r1, r0, n2, n1, n0, d1, d0) FLINT_MPN_UDIV_QR_3BY2_HW(q, r1, r0, n2, n1, n0, d1, d0)
    _FLINT_MPN_FDIV_3_NUMERATOR(r2, r1, r0, s, a1, a0, d2, d1, d0);
    _FLINT_MPN_FDIV_3_STEP(q2, r2, r1, r0, x, d2, d1, d0, _FDIV_3_DIV3BY2);
    _FLINT_MPN_FDIV_3_STEP(q1, r2, r1, r0, x, d2, d1, d0, _FDIV_3_DIV3BY2);
    _FLINT_MPN_FDIV_3_STEP(q0, r2, r1, r0, x, d2, d1, d0, _FDIV_3_DIV3BY2);
#undef _FDIV_3_DIV3BY2
    *q2p = q2;
    *q1p = q1;
    *q0p = q0;
    *sp = s;
    return (r2 | r1 | r0) != 0;
#else
    return _flint_mpn_fdiv_3_kernel_preinv(q2p, q1p, q0p, sp, a1, a0, d2, d1, d0, flint_mpn_preinv1(d2, d1));
#endif
}

/* The square root kernels below use double precision for the first
   approximation of a 64-bit limb and are not written for FLINT_BITS == 32. */
#if FLINT_BITS == 64

/*
    One limb. With c = 4A (t = 0) or 2A (t = 1), y = 1/sqrt(c) is computed
    in double precision: c is rounded once on conversion, then sqrt and
    division are correctly rounded, so |y / M - 1| < 1.25 * 2^-52; we take
    Y0 = floor(y 2^64) (clamped below 2^64), with relative error
    eta < 2^-51.6. With the residual eps = 1 - c Y0^2 / 2^128 evaluated
    exactly as E2 / 2^(190 + t), E2 = 2^(190 + t) - a Y0^2, we have
    M = y0 (1 - eps)^(-1/2) = y0 (1 + eps/2 + 3 eps^2 / 8 + 5 eps^3 / 16 + ...),
    |eps| < 2^-50.6. The second-order term y0 eps / 2 is evaluated at
    128-bit scale as Y0 E2 2^(-127 - t) (truncations of E2 / 2^64 and of
    the product: below 4 units of 2^-128 in either direction), and the
    third-order term (3/8) eps^2 y0, which is below 2^27 units, in double
    precision (error below 1 unit, plus 1 for the conversion); the
    remaining terms are below 2^-24 units. Thus |Y1 - M 2^128| < 6 and,
    with D the low limb of Y1, the high limb is the certified (inexact)
    truncation when 8 <= D <= 2^64 - 9.

    Returns 1 if certified (or if no certification is required), writing
    the truncated mantissa to *m.
*/
FLINT_FORCE_INLINE void
_flint_mpn_rsqrt_1_approx_128(ulong * Y1hi, ulong * Y1lo, ulong a, ulong t, int third_order)
{
    double c, y;
    ulong Y0, p1, p0, h1, l1, h0, l0, q2, q1, q0;
    ulong e1, e0, x1, x0, m1, m0, n1, n0, r2, r1, r0, d1, d0, y1, y0;
    int neg;

    c = (double) a * (t ? 0x1p-63 : 0x1p-62);
    y = 1.0 / sqrt(c);
    Y0 = (y >= 1.0) ? UWORD_MAX : (ulong) (y * 0x1p64);

    /* a Y0^2 as (q2, q1, q0) */
    umul_ppmm(p1, p0, Y0, Y0);
    umul_ppmm(h1, l1, a, p1);
    umul_ppmm(h0, l0, a, p0);
    q0 = l0;
    add_ssaaaa(q2, q1, h1, l1, 0, h0);

    /* (e1, e0) = floor(E2 / 2^64), signed */
    sub_ddmmss(e1, e0, UWORD(1) << (62 + t), 0, q2, q1);
    sub_ddmmss(e1, e0, e1, e0, 0, (q0 != 0));

    neg = ((slong) e1 < 0);
    if (neg)
        sub_ddmmss(x1, x0, 0, 0, e1, e0);
    else
    {
        x1 = e1;
        x0 = e0;
    }

    /* (r2, r1, r0) = Y0 |X| */
    umul_ppmm(m1, m0, Y0, x0);
    umul_ppmm(n1, n0, Y0, x1);
    r0 = m0;
    add_ssaaaa(r2, r1, n1, n0, 0, m1);

    /* shift right by 63 + t */
    if (t)
    {
        d1 = r2;
        d0 = r1;
    }
    else
    {
        d1 = (r2 << 1) | (r1 >> (FLINT_BITS - 1));
        d0 = (r1 << 1) | (r0 >> (FLINT_BITS - 1));
    }

    if (neg)
        sub_ddmmss(y1, y0, Y0, 0, d1, d0);
    else
        add_ssaaaa(y1, y0, Y0, 0, d1, d0);

    /* third-order term (3/8) eps^2 Y0 2^64 < 2^27, in double precision
       (only needed for certification: without it, |Y1 - M 2^128| < 2^26) */
    if (third_order)
    {
        double eps, c2;
        eps = ((double) x1 * 0x1p64 + (double) x0) * (t ? 0x1p-127 : 0x1p-126);
        c2 = 0.375 * eps * eps * ((double) Y0 * 0x1p64);
        add_ssaaaa(y1, y0, y1, y0, 0, (ulong) c2);
    }

    *Y1hi = y1;
    *Y1lo = y0;
}

FLINT_FORCE_INLINE int
_flint_mpn_rsqrt_1_approx(ulong * m, ulong a, ulong t, int certify)
{
    ulong y1, y0;

    _flint_mpn_rsqrt_1_approx_128(&y1, &y0, a, t, certify);
    *m = y1;

    if (!certify)
        return 1;

    return (y0 >= 8) && (y0 <= UWORD_MAX - 8);
}


/*
    Square roots at two and three fraction limbs (64-bit only).

    _flint_mpn_small_rsqrt_2: Y = (y1, y0) with |Y - M B^2| < 8 for
    M = (4A)^(-1/2) (t = 0) or (2A)^(-1/2) (t = 1), A in [1/2, 1) the
    fraction whose top two limbs are (a1, a0) (a1 >= B/2): the double
    precision 1/sqrt, refined at 128 bits with the residual of both limbs
    and the third-order term (see the comments below).

    _flint_mpn_small_sqrt_2: S = sqrt(A 2^-t) in [1/2, 1) at two fraction
    limbs as A M 2^(1 - t), within FLINT_MPN_SMALL_SQRT_ERR_2 units (the
    product below 3 units, 8 A from Y, one from the truncation of A; doubled
    for t = 0).

    _flint_mpn_small_sqrt_3: S at three fraction limbs for A = (a2, a1, a0),
    from S0 at two limbs and one Newton step S = S0 + R M with the residual
    R = A 2^-t - S0^2 computed exactly (|R| < 2^-122): the step leaves an
    error of order (24 B^-2)^2, M is within 2^-123 relatively, and the
    product R M is truncated, within FLINT_MPN_SMALL_SQRT_ERR_3 units.

    Both clamp S to [1/2, 1 - ulp].
*/
#define FLINT_MPN_SMALL_SQRT_ERR_2 24
#define FLINT_MPN_SMALL_SQRT_ERR_3 3

FLINT_FORCE_INLINE void
_flint_mpn_small_rsqrt_2(ulong * Yhi, ulong * Ylo, ulong a1, ulong a0, ulong t)
{
    double c, y;
    ulong Y0, p1, p0, w3, w2, w1, w0, e3, e2, e1, x2, x1, x0;
    ulong r3, r2, r1, r0, u1, u0, d1, d0, y1, y0, cy;
    int neg;

    c = (double) a1 * (t ? 0x1p-63 : 0x1p-62);
    y = 1.0 / sqrt(c);
    Y0 = (y >= 1.0) ? UWORD_MAX : (ulong) (y * 0x1p64);

    umul_ppmm(p1, p0, Y0, Y0);
    FLINT_MPN_MUL_2X2(w3, w2, w1, w0, a1, a0, p1, p0);

    /* (e3, e2, e1, e0) = 2^(254 + t) - a Y0^2 */
    sub_dddmmmsss(e3, e2, e1, UWORD(1) << (62 + t), 0, 0, w3, w2, w1);
    sub_dddmmmsss(e3, e2, e1, e3, e2, e1, 0, 0, (w0 != 0));

    /* |floor(E / 2^64)| as (x2, x1, x0) */
    neg = ((slong) e3 < 0);
    if (neg)
        sub_dddmmmsss(x2, x1, x0, 0, 0, 0, e3, e2, e1);
    else
    {
        x2 = e3;
        x1 = e2;
        x0 = e1;
    }
    FLINT_ASSERT(x2 < (UWORD(1) << 20));

    /* (r3, r2, r1, r0) = Y0 |X| */
    umul_ppmm(r1, r0, Y0, x0);
    umul_ppmm(u1, u0, Y0, x1);
    add_ssaaaa(r2, r1, u1, u0, 0, r1);
    r3 = 0;
    umul_ppmm(u1, u0, Y0, x2);
    add_ssaaaa(r3, r2, u1, u0, 0, r2);

    /* shift right by 127 + t bits */
    if (t)
    {
        d1 = r3;
        d0 = r2;
    }
    else
    {
        d1 = (r3 << 1) | (r2 >> (FLINT_BITS - 1));
        d0 = (r2 << 1) | (r1 >> (FLINT_BITS - 1));
    }

    /* third-order term (3/8) eps^2 Y0 2^64 < 2^27, in double precision */
    {
        double eps, c2;
        eps = ((double) x2 * 0x1p128 + (double) x1 * 0x1p64 + (double) x0) * (t ? 0x1p-191 : 0x1p-190);
        c2 = 0.375 * eps * eps * ((double) Y0 * 0x1p64);
        if (neg)
            sub_ddmmss(d1, d0, d1, d0, 0, (ulong) c2);
        else
            add_ssaaaa(d1, d0, d1, d0, 0, (ulong) c2);
    }

    if (neg)
    {
        sub_ddmmss(y1, y0, Y0, 0, d1, d0);
    }
    else
    {
        add_sssaaaaaa(cy, y1, y0, 0, Y0, 0, 0, d1, d0);
        if (cy)
        {
            y1 = UWORD_MAX;
            y0 = UWORD_MAX;
        }
    }

    *Yhi = y1;
    *Ylo = y0;
}

/* S = A M 2^(1 - t) from Y, clamped */
FLINT_FORCE_INLINE void
_flint_mpn_small_sqrt_from_rsqrt_2(ulong * s1, ulong * s0, ulong a1, ulong a0, ulong y1, ulong y0, ulong t)
{
    ulong p1, p0;

    {
        ulong l1, l0;
        FLINT_MPN_MUL_2X2(p1, p0, l1, l0, a1, a0, y1, y0);
        (void) l1; (void) l0;
    }
    if (t == 0)
    {
        if (p1 >> (FLINT_BITS - 1))
        {
            p1 = UWORD_MAX;
            p0 = UWORD_MAX;
        }
        else
        {
            p1 = (p1 << 1) | (p0 >> (FLINT_BITS - 1));
            p0 = p0 << 1;
        }
    }
    if (!(p1 >> (FLINT_BITS - 1)))
    {
        p1 = UWORD(1) << (FLINT_BITS - 1);
        p0 = 0;
    }
    *s1 = p1;
    *s0 = p0;
}

FLINT_FORCE_INLINE void
_flint_mpn_small_sqrt_2(ulong * s1, ulong * s0, ulong a1, ulong a0, ulong t)
{
    ulong y1, y0;
    _flint_mpn_small_rsqrt_2(&y1, &y0, a1, a0, t);
    _flint_mpn_small_sqrt_from_rsqrt_2(s1, s0, a1, a0, y1, y0, t);
}

FLINT_FORCE_INLINE void
_flint_mpn_small_sqrt_3(ulong * s2, ulong * s1, ulong * s0, ulong a2, ulong a1, ulong a0, ulong t)
{
    ulong y1, y0, t1, t0, q3, q2, q1, q0, A3, A2, A1, A0;
    ulong r2, r1, r0, p3, p2, p1, p0, h2, h1, h0, cy;
    int neg;

    _flint_mpn_small_rsqrt_2(&y1, &y0, a2, a1, t);
    _flint_mpn_small_sqrt_from_rsqrt_2(&t1, &t0, a2, a1, y1, y0, t);

    /* R B^4 = A 2^-t B^4 - S0^2 B^4, exactly (|R| B^4 < 2^134) */
    FLINT_MPN_SQR_2X2(q3, q2, q1, q0, t1, t0);
    if (t)
    {
        A3 = a2 >> 1;
        A2 = (a2 << (FLINT_BITS - 1)) | (a1 >> 1);
        A1 = (a1 << (FLINT_BITS - 1)) | (a0 >> 1);
        A0 = a0 << (FLINT_BITS - 1);
    }
    else
    {
        A3 = a2; A2 = a1; A1 = a0; A0 = 0;
    }
    sub_ddddmmmmssss(A3, A2, A1, A0, A3, A2, A1, A0, q3, q2, q1, q0);
    neg = ((slong) A3 < 0);
    if (neg)
        sub_ddddmmmmssss(A3, A2, A1, A0, UWORD(0), UWORD(0), UWORD(0), UWORD(0), A3, A2, A1, A0);
    r2 = A2; r1 = A1; r0 = A0;

    /* floor(|R| B^4 Y / B^3) at two limbs: delta B^3 = |R| M B^3 */
    FLINT_MPN_MUL_2X2(p3, p2, p1, p0, r1, r0, y1, y0);
    FLINT_MPN_MUL_2X1(h2, h1, h0, y1, y0, r2);
    add_sssaaaaaa(cy, p3, p2, UWORD(0), p3, p2, UWORD(0), h1, h0);
    h2 += cy;
    (void) p1; (void) p0;

    /* S = S0 B + delta (h2, p3) */
    if (neg)
    {
        sub_dddmmmsss(t1, t0, q0, t1, t0, UWORD(0), UWORD(0), h2, p3);
    }
    else
    {
        q0 = 0;
        add_ssssaaaaaaaa(cy, t1, t0, q0, UWORD(0), t1, t0, q0, UWORD(0), UWORD(0), h2, p3);
        if (cy)
        {
            t1 = t0 = q0 = UWORD_MAX;
        }
    }

    if (!(t1 >> (FLINT_BITS - 1)))
    {
        t1 = UWORD(1) << (FLINT_BITS - 1);
        t0 = q0 = 0;
    }

    *s2 = t1;
    *s1 = t0;
    *s0 = q0;
}

/*
    _flint_mpn_small_rsqrt_3: Z = (z2, z1, z0) with |Z - M B^3| < 2 for the
    exact two-limb fraction A = (a1, a0) (a1 >= B/2), M as above, by one
    Newton step from Y = (y1, y0) of _flint_mpn_small_rsqrt_2: with y = Y B^-2
    and c = 2^(2 - t) A, z = y + y eps / 2 where eps = 1 - c y^2.

    Since |y - M| < 8 B^-2 and M > 1/2, y = M (1 + delta) with
    |delta| < 2^-124, and the exact step leaves M (1.5 delta^2 + ...),
    below 2^-55 units of B^-3. The residual E = B^6 eps = B^6 - 2^(2 - t) P,
    P = A Y^2 B^? exactly (six limbs), has |E| < 2^262; it is used down to
    the limb of weight B^2 only (an error below B^2, i.e. half a unit in the
    correction Y E / (2 B^5)), and the correction is truncated (below one
    unit). Results reaching B^3 (only for t = 1, M close to 1) are clamped
    to B^3 - 1, which keeps the bound since M B^3 < B^3.
*/
#define FLINT_MPN_SMALL_RSQRT_ERR_3 2

/*
    s = floor(sqrt(N)) exactly for N = (n1, n0) with n1 >= B/4 (so that
    s >= B/2), returning whether N is not a perfect square; the _rem
    variant also gives the remainder N - s^2 = (rem1, rem0) <= 2 s.

    The double precision root y = sqrt(N) (1 + epsilon), |epsilon| < 2^-52
    (N converted with relative error below 2^-52, halved by the square root,
    which is correctly rounded), gives s0 = floor(y) within
    2^12 + 1 of sqrt(N); the exact residual r = N - s0^2, |r| < 2^78, in
    double precision divided by 2 y gives the Newton correction within
    2^-30 (relative errors of a few 2^-53 on a quantity below 2^13, and the
    Newton error (sqrt(N) - s0)^2 / (2 s0) < 2^-38), so that
    s0 + floor(r / (2 y)) is floor(sqrt(N)) or one off, which the exact
    residual of the result corrects. No overflow: sqrt(N) < B - 1/2.
*/
FLINT_FORCE_INLINE int
_flint_mpn_isqrt_2_rem_inv(ulong * s, ulong * rem1, ulong * rem0, double * inv2y, ulong n1, ulong n0)
{
    double y, rd, inv;
    ulong s0, p1, p0, r1, r0;

    y = sqrt((double) n1 * 0x1p64 + (double) n0);
    s0 = (y >= 0x1p64) ? UWORD_MAX : (ulong) y;

    /* 1 / (2 y), also returned: a remainder times it is the fraction
       sqrt(N) - s to double precision (y and s differ by less than one
       part in 2^63) */
    inv = 0.5 / y;
    *inv2y = inv;

    umul_ppmm(p1, p0, s0, s0);
    sub_ddmmss(r1, r0, n1, n0, p1, p0);
    rd = (double) (slong) r1 * 0x1p64 + (double) r0;
    s0 += (slong) floor(rd * inv);

    /* exact residual of s0, then at most one step either way */
    umul_ppmm(p1, p0, s0, s0);
    sub_ddmmss(r1, r0, n1, n0, p1, p0);
    if ((slong) r1 < 0)
    {
        /* (s0 - 1)^2 = s0^2 - 2 s0 + 1 */
        add_ssaaaa(r1, r0, r1, r0, UWORD(0), s0);
        add_ssaaaa(r1, r0, r1, r0, UWORD(0), s0 - 1);
        s0--;
    }
    else
    {
        /* (s0 + 1)^2 = s0^2 + 2 s0 + 1, with 2 s0 + 1 = (t1, t0) */
        ulong t1 = s0 >> (FLINT_BITS - 1), t0 = (s0 << 1) | 1;

        if (r1 > t1 || (r1 == t1 && r0 >= t0))
        {
            sub_ddmmss(r1, r0, r1, r0, t1, t0);
            s0++;
        }
    }
    *s = s0;
    *rem1 = r1;
    *rem0 = r0;
    return (r1 | r0) != 0;
}

FLINT_FORCE_INLINE int
_flint_mpn_isqrt_2_rem(ulong * s, ulong * rem1, ulong * rem0, ulong n1, ulong n0)
{
    double inv;
    return _flint_mpn_isqrt_2_rem_inv(s, rem1, rem0, &inv, n1, n0);
}

/* the same without the remainder N - s^2 (which is below 2 s + 1) */
FLINT_FORCE_INLINE int
_flint_mpn_isqrt_2(ulong * s, ulong n1, ulong n0)
{
    ulong r1, r0;
    return _flint_mpn_isqrt_2_rem(s, &r1, &r0, n1, n0);
}

FLINT_FORCE_INLINE void
_flint_mpn_small_rsqrt_3(ulong * z2, ulong * z1, ulong * z0, ulong a1, ulong a0, ulong t)
{
    ulong y1, y0, q3, q2, q1, q0, l3, l2, l1, l0, h3, h2, h1, h0;
    ulong p5, p4, p3, p2, p1, g6, g5, g4, g3, g2, x2, x1, x0;
    ulong m3, m2, m1, m0, k2, k1, k0, w4, w3, c1, c0, cy;
    unsigned int s = 2 - t;
    int neg;

    _flint_mpn_small_rsqrt_2(&y1, &y0, a1, a0, t);

    /* P = A Y^2 (B^6 scale), exactly; only limbs 1 to 5 are needed */
    FLINT_MPN_SQR_2X2(q3, q2, q1, q0, y1, y0);
    FLINT_MPN_MUL_2X2(l3, l2, l1, l0, a1, a0, q1, q0);
    FLINT_MPN_MUL_2X2(h3, h2, h1, h0, a1, a0, q3, q2);
    (void) l0;
    p1 = l1;
    add_ssssaaaaaaaa(p5, p4, p3, p2, h3, h2, h1, h0, UWORD(0), UWORD(0), l3, l2);

    /* G = 2^(2 - t) P, limbs 2 to 6 */
    g6 = p5 >> (FLINT_BITS - s);
    g5 = (p5 << s) | (p4 >> (FLINT_BITS - s));
    g4 = (p4 << s) | (p3 >> (FLINT_BITS - s));
    g3 = (p3 << s) | (p2 >> (FLINT_BITS - s));
    g2 = (p2 << s) | (p1 >> (FLINT_BITS - s));
    (void) g5;    /* only checked by the assertions below */

    /* |E| / B^2 as (x2, x1, x0), E = B^6 - G with the limbs of G below
       B^2 dropped */
    if (g6 != 0)
    {
        neg = 1;
        FLINT_ASSERT(g6 == 1 && g5 == 0);
        x2 = g4; x1 = g3; x0 = g2;
    }
    else
    {
        neg = 0;
        FLINT_ASSERT(g5 == UWORD_MAX);
        sub_dddmmmsss(x2, x1, x0, UWORD(0), UWORD(0), UWORD(0), g4, g3, g2);
    }
    FLINT_ASSERT(x2 < (UWORD(1) << 8));

    /* (c1, c0) = floor(Y |X| / (2 B^3)) */
    FLINT_MPN_MUL_2X2(m3, m2, m1, m0, y1, y0, x1, x0);
    FLINT_MPN_MUL_2X1(k2, k1, k0, y1, y0, x2);
    (void) m0; (void) m1;
    add_sssaaaaaa(w4, w3, cy, k2, k1, k0, UWORD(0), m3, m2);
    (void) cy;
    c1 = w4 >> 1;
    c0 = (w3 >> 1) | (w4 << (FLINT_BITS - 1));

    /* Z = Y B -+ (c1, c0) */
    if (neg)
    {
        sub_dddmmmsss(y1, y0, k0, y1, y0, UWORD(0), UWORD(0), c1, c0);
    }
    else
    {
        add_ssssaaaaaaaa(cy, y1, y0, k0, UWORD(0), y1, y0, UWORD(0), UWORD(0), UWORD(0), c1, c0);
        if (cy)
            y1 = y0 = k0 = UWORD_MAX;
    }

    *z2 = y1;
    *z1 = y0;
    *z0 = k0;
}
#endif /* FLINT_BITS == 64 */

#endif
