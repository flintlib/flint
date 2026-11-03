/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Internal definitions for the elementary functions. */

#ifndef NFLOAT_ELEM_H
#define NFLOAT_ELEM_H

#include <math.h>
#include "longlong.h"
#include "mpn_extras.h"
#include "mp_real.h"
#include "nfloat.h"
#include "impl.h"
#include "../mp_real/impl.h"
#include "../mpn_extras/fp_small.h"

/*
    Conventions.

    The functions evaluate fixed-point approximations with the kernels of
    the mp_real module and an error bound in ulps of the fixed-point
    result, which _nfloat_set_mpn_err turns into an nfloat: truncated in
    the default rounding mode, and in the floor and ceiling modes the
    endpoint of the error interval in the rounding direction, rounded
    outward, so that the result is a valid bound.

    The target is a relative accuracy of 2^-p where p is the context's
    function precision (NFLOAT_CTX_FUNC_PREC, by default the full
    precision); the kernels run at enough limbs for p plus
    NFLOAT_ELEM_GUARD_BITS bits, which covers their error bounds (below
    2^11 ulps, under 2^-1 ulp at p bits) and the error of the argument
    reduction.
*/

#define NFLOAT_ELEM_GUARD_BITS 12

/* Limbs of temporaries: target precision plus guard bits plus as many
   bits again for arguments close to zeros of the function. */
#define NFLOAT_ELEM_MAX_LIMBS (2 * NFLOAT_MAX_LIMBS + 8)

/* The largest left shift in normalizing a fixed-point value of small
   magnitude to a floating-point mantissa, which scales up its error count:
   that must stay below 2^(FLINT_BITS - 4) (the bound of the _nfloat_rf_t
   operations) for kernel errors below 2^11 ulps and a few operations, so
   40 bits on 64-bit machines and 8 on 32-bit machines; values with more
   leading zero bits are computed at more limbs. */
#define NFLOAT_ELEM_MAX_NORM_SHIFT (FLINT_BITS - 24)

/* fixed-point limbs for an absolute accuracy of 2^-(p + extra) */
FLINT_FORCE_INLINE slong
_nfloat_elem_limbs(gr_ctx_t ctx, slong extra)
{
    slong n = (NFLOAT_CTX_FUNC_PREC(ctx) + NFLOAT_ELEM_GUARD_BITS + extra
        + FLINT_BITS - 1) / FLINT_BITS;
#if FLINT_BITS == 32
    /* some kernels require two limbs on 32-bit machines */
    n = FLINT_MAX(n, 2);
#endif
    return n;
}

/* Sets (v, vn) to floor(|x| B^fn) mod B^vn for the nfloat
   |x| = (d, n) 2^(e - FLINT_BITS n) (the caller ensures |x| < B^(vn - fn)
   when the value matters in full); returns 1 if nonzero bits of |x| B^fn
   were dropped below the unit, otherwise 0. */
FLINT_FORCE_INLINE int
_nfloat_get_fixed(nn_ptr v, slong vn, slong fn, nn_srcptr d, slong n, slong e)
{
    slong s = e + FLINT_BITS * (fn - n);
    slong ls, j;
    unsigned int bs;
    int inexact;

    if (s >= 0)
    {
        ls = s / FLINT_BITS;
        bs = s % FLINT_BITS;

        /* (explicit helpers: the compiler would turn plain loops into
           calls to memset and memcpy, expensive for a few limbs) */
        _nfloat_zero_limbs(v, FLINT_MIN(ls, vn));

        if (ls < vn)
        {
            slong m = FLINT_MIN(n, vn - ls);

            if (bs == 0)
            {
                _nfloat_copy_limbs(v + ls, d, m);
                _nfloat_zero_limbs(v + ls + m, vn - ls - m);
            }
            else
            {
                v[ls] = d[0] << bs;
                for (j = 1; j < m; j++)
                    v[ls + j] = (d[j] << bs) | (d[j - 1] >> (FLINT_BITS - bs));
                if (ls + m < vn)
                {
                    v[ls + m] = d[m - 1] >> (FLINT_BITS - bs);
                    _nfloat_zero_limbs(v + ls + m + 1, vn - ls - m - 1);
                }
            }
        }

        return 0;
    }
    else
    {
        s = -s;
        ls = s / FLINT_BITS;
        bs = s % FLINT_BITS;

        if (ls >= n)
        {
            _nfloat_zero_limbs(v, vn);
            return 1;
        }

        inexact = 0;
        for (j = 0; j < ls; j++)
            inexact |= (d[j] != 0);

        if (bs == 0)
        {
            for (j = 0; j < vn; j++)
                v[j] = (ls + j < n) ? d[ls + j] : 0;
        }
        else
        {
            inexact |= ((d[ls] << (FLINT_BITS - bs)) != 0);
            for (j = 0; j < vn; j++)
            {
                ulong lo = (ls + j < n) ? d[ls + j] : 0;
                ulong hi = (ls + j + 1 < n) ? d[ls + j + 1] : 0;
                v[j] = (lo >> bs) | (hi << (FLINT_BITS - bs));
            }
        }

        return inexact;
    }
}

/* Sets res to (-1)^sgnbit V where V is a real number with
   (Y - elo) 2^e <= V <= (Y + ehi) 2^e for the integer Y = (y, len) > 0
   (y must have room for len + 1 limbs; it is destroyed). */
int _nfloat_set_mpn_err(nfloat_ptr res, nn_ptr y, slong len, slong e, ulong elo, ulong ehi, int sgnbit, gr_ctx_t ctx);

#if FLINT_BITS == 64
/* The finishing step of the in-register fast paths at one and two limbs
   (64-bit only): res = (-1)^sgnbit Y 2^(k - FLINT_BITS (len - 1)) for
   Y = (y2, y1, y0) resp. (y3, y2, y1, y0) with the integral limb y2 resp.
   y3 below 2^(FLINT_BITS - 1), the top two limbs not both zero, known
   within err units: truncated to the context's limbs in the default
   rounding mode, through _nfloat_set_mpn_err otherwise. */
FLINT_FORCE_INLINE int
_nfloat_small_set_1(nfloat_ptr res, ulong y2, ulong y1, ulong y0, slong k,
    ulong err, int sgnbit, gr_ctx_t ctx)
{
    unsigned int c;

    if (FLINT_UNLIKELY(NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx)))
    {
        ulong t[4];
        t[0] = y0; t[1] = y1; t[2] = y2;
        return _nfloat_set_mpn_err(res, t, 3, k - 2 * FLINT_BITS, err, err, sgnbit, ctx);
    }

    if (y2 != 0)
    {
        c = flint_clz(y2);
        NFLOAT_D(res)[0] = (y2 << c) | (y1 >> (FLINT_BITS - c));
        NFLOAT_EXP(res) = k + FLINT_BITS - c;
    }
    else
    {
        c = flint_clz(y1);
        NFLOAT_D(res)[0] = (c == 0) ? y1 : (y1 << c) | (y0 >> (FLINT_BITS - c));
        NFLOAT_EXP(res) = k - (slong) c;
    }
    NFLOAT_SGNBIT(res) = sgnbit;
    NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
    return GR_SUCCESS;
}

FLINT_FORCE_INLINE int
_nfloat_small_set_2(nfloat_ptr res, ulong y3, ulong y2, ulong y1, ulong y0, slong k,
    ulong err, int sgnbit, gr_ctx_t ctx)
{
    unsigned int c;

    if (FLINT_UNLIKELY(NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx)))
    {
        ulong t[5];
        t[0] = y0; t[1] = y1; t[2] = y2; t[3] = y3;
        return _nfloat_set_mpn_err(res, t, 4, k - 3 * FLINT_BITS, err, err, sgnbit, ctx);
    }

    if (y3 != 0)
    {
        c = flint_clz(y3);
        NFLOAT_D(res)[1] = (y3 << c) | (y2 >> (FLINT_BITS - c));
        NFLOAT_D(res)[0] = (y2 << c) | (y1 >> (FLINT_BITS - c));
        NFLOAT_EXP(res) = k + FLINT_BITS - c;
    }
    else
    {
        c = flint_clz(y2);
        if (c == 0)
        {
            NFLOAT_D(res)[1] = y2;
            NFLOAT_D(res)[0] = y1;
        }
        else
        {
            NFLOAT_D(res)[1] = (y2 << c) | (y1 >> (FLINT_BITS - c));
            NFLOAT_D(res)[0] = (y1 << c) | (y0 >> (FLINT_BITS - c));
        }
        NFLOAT_EXP(res) = k - (slong) c;
    }
    NFLOAT_SGNBIT(res) = sgnbit;
    NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
    return GR_SUCCESS;
}

/* the return value of a fast path that does not apply to its arguments */
#define NFLOAT_FAST_DECLINE (-1)

/* the cotangent and the hyperbolic tangent and cotangent by the fast
   paths (sin_cos.c, hyperbolic.c) for one and two limbs, or
   NFLOAT_FAST_DECLINE */
#define NFLOAT_TRIG_COT 1
#define NFLOAT_TRIG_SEC 2
#define NFLOAT_TRIG_CSC 3
#define NFLOAT_TRIG_SINC 4
int _nfloat_trig_recip_fast(nfloat_ptr res, nfloat_srcptr x, int pi, int which, gr_ctx_t ctx);

/* inverse functions at one and two limbs (elem_arc_fast.c), for normal x
   in the domain (|x| != 1 where that is special) */
#define NFLOAT_ARC_ASIN 0
#define NFLOAT_ARC_ACOS 1
#define NFLOAT_ARC_ATAN 2
#define NFLOAT_ARC_ACOT 3
#define NFLOAT_ARC_ASEC 4
#define NFLOAT_ARC_ACSC 5
#define NFLOAT_ARC_ASINH 6
#define NFLOAT_ARC_ACOSH 7
#define NFLOAT_ARC_ATANH 8
#define NFLOAT_ARC_ACOTH 9
#define NFLOAT_ARC_ASECH 10
#define NFLOAT_ARC_ACSCH 11
int _nfloat_arc_fast(nfloat_ptr res, nfloat_srcptr x, int f, int div_pi, gr_ctx_t ctx);
int _nfloat_tanh_coth_fast(nfloat_ptr res, nfloat_srcptr x, int coth, gr_ctx_t ctx);

/* Fixed-point helpers of the fast paths (64-bit). */

/* |x| = d 2^(e - FLINT_BITS) as (xi, x1, x0) with the binary point after
   the integral limb xi, exactly, for -FLINT_BITS <= e < FLINT_BITS */
FLINT_FORCE_INLINE void
_nfloat_small_fixed_1(ulong * xi, ulong * x1, ulong * x0, ulong d, slong e)
{
    if (e > 0)
    {
        *xi = d >> (FLINT_BITS - e);
        *x1 = d << e;
        *x0 = 0;
    }
    else if (e == 0)
    {
        *xi = 0; *x1 = d; *x0 = 0;
    }
    else if (e > -FLINT_BITS)
    {
        *xi = 0;
        *x1 = d >> (-e);
        *x0 = d << (FLINT_BITS + e);
    }
    else
    {
        *xi = 0; *x1 = 0; *x0 = d;
    }
}

/* |x| = (d1, d0) 2^(e - 2 FLINT_BITS) as (xi, x2, x1, x0), binary point
   after xi, for -2 FLINT_BITS <= e < FLINT_BITS: exact for e >=
   -FLINT_BITS, truncated below B^-3 otherwise */
FLINT_FORCE_INLINE void
_nfloat_small_fixed_2(ulong * xi, ulong * x2, ulong * x1, ulong * x0, ulong d1, ulong d0, slong e)
{
    unsigned int s;

    *xi = *x2 = *x1 = *x0 = 0;
    if (e > 0)
    {
        *xi = d1 >> (FLINT_BITS - e);
        *x2 = (d1 << e) | (d0 >> (FLINT_BITS - e));
        *x1 = d0 << e;
    }
    else if (e == 0)
    {
        *x2 = d1; *x1 = d0;
    }
    else if (e > -FLINT_BITS)
    {
        s = -e;
        *x2 = d1 >> s;
        *x1 = (d1 << (FLINT_BITS - s)) | (d0 >> s);
        *x0 = d0 << (FLINT_BITS - s);
    }
    else if (e == -FLINT_BITS)
    {
        *x1 = d1; *x0 = d0;
    }
    else if (e > -2 * FLINT_BITS)
    {
        s = -e - FLINT_BITS;
        *x1 = d1 >> s;
        *x0 = (d1 << (FLINT_BITS - s)) | (d0 >> s);
    }
    else
        *x0 = d1;
}

/* The reduction of exp: for X = (xi, x2, x1, x0) >= 0 with three fraction
   limbs and xi < 2^62, returns q and sets t = X - q L in [0, L) (three
   fraction limbs), L the floor of log 2 at three limbs: q from
   _mp_real_exp_quotient (exact or one low, corrected once); t is within
   (q + 1) 2^-192 of X - q log 2. */
FLINT_FORCE_INLINE ulong
_nfloat_small_exp_reduce_3(ulong * t2, ulong * t1, ulong * t0,
    ulong xi, ulong x2, ulong x1, ulong x0)
{
    nn_srcptr L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, 3);
    ulong q, p3, p2, p1, p0, h, l, u3, u2, u1, u0;

    q = _mp_real_exp_quotient(xi, x2);
    umul_ppmm(p1, p0, q, L[0]);
    umul_ppmm(h, l, q, L[1]);
    add_ssaaaa(p2, p1, h, l, UWORD(0), p1);
    umul_ppmm(h, l, q, L[2]);
    add_ssaaaa(p3, p2, h, l, UWORD(0), p2);
    sub_ddddmmmmssss(u3, u2, u1, u0, xi, x2, x1, x0, p3, p2, p1, p0);

    if (u3 != 0 || u2 > L[2] || (u2 == L[2] && (u1 > L[1] || (u1 == L[1] && u0 >= L[0]))))
    {
        sub_ddddmmmmssss(u3, u2, u1, u0, u3, u2, u1, u0, UWORD(0), L[2], L[1], L[0]);
        q++;
    }

    *t2 = u2; *t1 = u1; *t0 = u0;
    return q;
}

/* the same with four fraction limbs and L at four limbs (within (q + 1)
   2^-256) */
FLINT_FORCE_INLINE ulong
_nfloat_small_exp_reduce_4(ulong * t3, ulong * t2, ulong * t1, ulong * t0,
    ulong xi, ulong x3, ulong x2, ulong x1, ulong x0)
{
    nn_srcptr L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, 4);
    ulong q, p4, p3, p2, p1, p0, h, l, u4, u3, u2, u1, u0;

    q = _mp_real_exp_quotient(xi, x3);
    umul_ppmm(p1, p0, q, L[0]);
    umul_ppmm(h, l, q, L[1]);
    add_ssaaaa(p2, p1, h, l, UWORD(0), p1);
    umul_ppmm(h, l, q, L[2]);
    add_ssaaaa(p3, p2, h, l, UWORD(0), p2);
    umul_ppmm(h, l, q, L[3]);
    add_ssaaaa(p4, p3, h, l, UWORD(0), p3);
    sub_dddddmmmmmsssss(u4, u3, u2, u1, u0, xi, x3, x2, x1, x0, p4, p3, p2, p1, p0);

    if (u4 != 0 || u3 > L[3] || (u3 == L[3] && (u2 > L[2] || (u2 == L[2]
            && (u1 > L[1] || (u1 == L[1] && u0 >= L[0]))))))
    {
        sub_dddddmmmmmsssss(u4, u3, u2, u1, u0, u4, u3, u2, u1, u0, UWORD(0), L[3], L[2], L[1], L[0]);
        q++;
    }

    *t3 = u3; *t2 = u2; *t1 = u1; *t0 = u0;
    return q;
}

/* exp(x) at two resp. three limbs from the reduction above and the
   table-driven kernel (exp.c explains the error bounds) */
#define NFLOAT_SMALL_EXP_ERR_2 (MP_REAL_SMALL_EXP_ERR_2 + 3)
#define NFLOAT_SMALL_EXP_ERR_3 (MP_REAL_SMALL_EXP_ERR_3 + 3)

/* exp(x) = Y 2^(k - 2 FLINT_BITS), Y = (y2, y1, y0) in [B^2, 2 B^2), within
   NFLOAT_SMALL_EXP_ERR_2 units, for the one-limb mantissa d */
FLINT_FORCE_INLINE slong
_nfloat_small_exp_1(ulong * y2, ulong * y1, ulong * y0, ulong d, slong e, int sgnbit)
{
    ulong xi, x1, x0, q, t2, t1, t0;
    slong k;

    _nfloat_small_fixed_1(&xi, &x1, &x0, d, e);
    q = _nfloat_small_exp_reduce_3(&t2, &t1, &t0, xi, x1, x0, 0);

    if (sgnbit)
    {
        nn_srcptr L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, 3);
        sub_dddmmmsss(t2, t1, t0, L[2], L[1], L[0], t2, t1, t0);
        k = -(slong) (q + 1);
    }
    else
        k = (slong) q;

    _mp_real_small_exp_2(y2, y1, y0, t2, t1);
    return k;
}

/* the same at three limbs: Y = (y3, y2, y1, y0), exp(x) = Y 2^(k - 3 FLINT_BITS) */
FLINT_FORCE_INLINE slong
_nfloat_small_exp_2(ulong * y3, ulong * y2, ulong * y1, ulong * y0, ulong d1, ulong d0, slong e, int sgnbit)
{
    ulong xi, x2, x1, x0, q, t3, t2, t1, t0;
    slong k;

    _nfloat_small_fixed_2(&xi, &x2, &x1, &x0, d1, d0, e);
    q = _nfloat_small_exp_reduce_4(&t3, &t2, &t1, &t0, xi, x2, x1, x0, 0);

    if (sgnbit)
    {
        nn_srcptr L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, 4);
        sub_ddddmmmmssss(t3, t2, t1, t0, L[3], L[2], L[1], L[0], t3, t2, t1, t0);
        k = -(slong) (q + 1);
    }
    else
        k = (slong) q;

    _mp_real_small_exp_3(y3, y2, y1, y0, t3, t2, t1);
    return k;
}

/* |log u| for u = M 2^E, M = (m1, m0) a normalized two-limb fraction and
   |E - 1| < 2^62: c L +- F with c = E - 1, F = log1p(2M - 1) by the
   table-driven kernel at two limbs (2M - 1 exact) and L the floor of
   log 2 at three fraction limbs (c times its error below a quarter ulp of
   B^-2), as (r2, r1, r0) with the integral limb r2, the lowest limb of
   the sum dropped: within NFLOAT_SMALL_LOG_ERR_2 units of B^-2. Returns
   1 if log u < 0 (c < 0). */
#define NFLOAT_SMALL_LOG_ERR_2 (MP_REAL_SMALL_LOG_ERR_2 + 2)
#define NFLOAT_SMALL_LOG_ERR_3 (MP_REAL_SMALL_LOG_ERR_3 + 2)

FLINT_FORCE_INLINE int
_nfloat_small_log_2(ulong * r2, ulong * r1, ulong * r0, ulong m1, ulong m0, slong E)
{
    nn_srcptr L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, 3);
    ulong f1, f0, cc, p3, p2, p1, p0, h, l;
    slong c = E - 1;

    _mp_real_small_log1p_2(&f1, &f0, (m1 << 1) | (m0 >> (FLINT_BITS - 1)), m0 << 1);

    cc = (c < 0) ? -(ulong) c : (ulong) c;
    umul_ppmm(p1, p0, cc, L[0]);
    umul_ppmm(h, l, cc, L[1]);
    add_ssaaaa(p2, p1, h, l, UWORD(0), p1);
    umul_ppmm(h, l, cc, L[2]);
    add_ssaaaa(p3, p2, h, l, UWORD(0), p2);
    (void) p0;

    if (c >= 0)
        add_sssaaaaaa(p3, p2, p1, p3, p2, p1, UWORD(0), f1, f0);
    else
        sub_dddmmmsss(p3, p2, p1, p3, p2, p1, UWORD(0), f1, f0);

    *r2 = p3; *r1 = p2; *r0 = p1;
    return c < 0;
}

/* the same at three limbs: M = (m2, m1, m0), the result (r3, r2, r1, r0)
   with three fraction limbs, within NFLOAT_SMALL_LOG_ERR_3 units of B^-3 */
FLINT_FORCE_INLINE int
_nfloat_small_log_3(ulong * r3, ulong * r2, ulong * r1, ulong * r0,
    ulong m2, ulong m1, ulong m0, slong E)
{
    nn_srcptr L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, 4);
    ulong f2, f1, f0, cc, p4, p3, p2, p1, p0, h, l;
    slong c = E - 1;

    _mp_real_small_log1p_3(&f2, &f1, &f0, (m2 << 1) | (m1 >> (FLINT_BITS - 1)),
        (m1 << 1) | (m0 >> (FLINT_BITS - 1)), m0 << 1);

    cc = (c < 0) ? -(ulong) c : (ulong) c;
    umul_ppmm(p1, p0, cc, L[0]);
    umul_ppmm(h, l, cc, L[1]);
    add_ssaaaa(p2, p1, h, l, UWORD(0), p1);
    umul_ppmm(h, l, cc, L[2]);
    add_ssaaaa(p3, p2, h, l, UWORD(0), p2);
    umul_ppmm(h, l, cc, L[3]);
    add_ssaaaa(p4, p3, h, l, UWORD(0), p3);
    (void) p0;

    if (c >= 0)
        add_ssssaaaaaaaa(p4, p3, p2, p1, p4, p3, p2, p1, UWORD(0), f2, f1, f0);
    else
        sub_ddddmmmmssss(p4, p3, p2, p1, p4, p3, p2, p1, UWORD(0), f2, f1, f0);

    *r3 = p4; *r2 = p3; *r1 = p2; *r0 = p1;
    return c < 0;
}

/* (y2, y1, y0) >> s for 0 <= s, truncated (zero for s >= 3 FLINT_BITS) */
FLINT_FORCE_INLINE void
_nfloat_small_rshift_3(ulong * r2, ulong * r1, ulong * r0, ulong y2, ulong y1, ulong y0, ulong s)
{
    unsigned int b;

    if (s >= 2 * FLINT_BITS)
    {
        *r2 = 0; *r1 = 0;
        *r0 = (s >= 3 * FLINT_BITS) ? 0 : y2 >> (s - 2 * FLINT_BITS);
        return;
    }
    if (s >= FLINT_BITS)
    {
        y0 = y1; y1 = y2; y2 = 0;
        s -= FLINT_BITS;
    }
    b = s;
    if (b != 0)
    {
        y0 = (y0 >> b) | (y1 << (FLINT_BITS - b));
        y1 = (y1 >> b) | (y2 << (FLINT_BITS - b));
        y2 = y2 >> b;
    }
    *r2 = y2; *r1 = y1; *r0 = y0;
}

/* (y3, y2, y1, y0) >> s, truncated */
FLINT_FORCE_INLINE void
_nfloat_small_rshift_4(ulong * r3, ulong * r2, ulong * r1, ulong * r0,
    ulong y3, ulong y2, ulong y1, ulong y0, ulong s)
{
    unsigned int b;

    if (s >= 4 * FLINT_BITS)
    {
        *r3 = *r2 = *r1 = *r0 = 0;
        return;
    }
    while (s >= FLINT_BITS)
    {
        y0 = y1; y1 = y2; y2 = y3; y3 = 0;
        s -= FLINT_BITS;
    }
    b = s;
    if (b != 0)
    {
        y0 = (y0 >> b) | (y1 << (FLINT_BITS - b));
        y1 = (y1 >> b) | (y2 << (FLINT_BITS - b));
        y2 = (y2 >> b) | (y3 << (FLINT_BITS - b));
        y3 = y3 >> b;
    }
    *r3 = y3; *r2 = y2; *r1 = y1; *r0 = y0;
}

/* the square root kernels now live in mpn_extras/fp_small.h */
#define _nfloat_rsqrt_1_approx_128 _flint_mpn_rsqrt_1_approx_128
#define _nfloat_rsqrt_1_approx _flint_mpn_rsqrt_1_approx
#define _nfloat_small_rsqrt_2 _flint_mpn_small_rsqrt_2
#define _nfloat_small_sqrt_2 _flint_mpn_small_sqrt_2
#define _nfloat_small_sqrt_3 _flint_mpn_small_sqrt_3
#define _nfloat_small_rsqrt_3 _flint_mpn_small_rsqrt_3
#define _nfloat_isqrt_2 _flint_mpn_isqrt_2
#define NFLOAT_SMALL_SQRT_ERR_2 FLINT_MPN_SMALL_SQRT_ERR_2
#define NFLOAT_SMALL_SQRT_ERR_3 FLINT_MPN_SMALL_SQRT_ERR_3
#define NFLOAT_SMALL_RSQRT_ERR_3 FLINT_MPN_SMALL_RSQRT_ERR_3

#endif

/* Sets res to (-1)^sgnbit (1 + delta) where 0 <= delta < 2^e (dsgn = 0)
   resp. -2^e < delta <= 0 (dsgn = 1), e <= -FLINT_BITS * nlimbs: that is,
   +-1 in the default rounding mode and a valid bound otherwise. */
int _nfloat_one_plus_tiny(nfloat_ptr res, int sgnbit, int dsgn, slong e, gr_ctx_t ctx);

/* Sets res to x (1 + delta) with -2^e <= delta <= 0 (side = -1) resp.
   0 <= delta <= 2^e (side = 1), e < -FLINT_BITS * nlimbs (a tiny
   relative perturbation): x in the default rounding mode, a valid bound
   otherwise. */
int _nfloat_set_tiny_rel_err(nfloat_ptr res, nfloat_srcptr x, int side, gr_ctx_t ctx);

/* arguments beyond 2^NFLOAT_EXP_MAX_ARG_EXP in absolute value overflow
   or underflow exp: exp(2^(FLINT_BITS - 3)) = 2^(2^(FLINT_BITS - 3) / log 2)
   is beyond the exponent range 2^(FLINT_BITS - 3) = NFLOAT_MAX_EXP + 1 */
#define NFLOAT_EXP_MAX_ARG_EXP (FLINT_BITS - 3)

/* exp(x) for the normal nfloat x = (-1)^sgnbit (d, n) 2^(e - FLINT_BITS n),
   e <= NFLOAT_EXP_MAX_ARG_EXP, with nk kernel limbs: sets (y, nk + 1)
   (y needs nk + 2 limbs) and *err, and returns k such that
   exp(x) = Y 2^(k - FLINT_BITS nk) within *err units, Y = (y, nk + 1). */
slong _nfloat_exp_mpn(nn_ptr y, ulong * err, nn_srcptr d, slong n, slong e, int sgnbit, slong nk);

/* A fixed-point result (d, len) 2^e within err units, with a sign, as
   produced by internal functions for use in compositions. */
typedef struct
{
    ulong d[NFLOAT_ELEM_MAX_LIMBS + 4];
    slong len;
    slong e;
    ulong err;
    int sgnbit;
}
_nfloat_fix_struct;

typedef _nfloat_fix_struct _nfloat_fix_t[1];

/* log(u) for u = (d, dn) 2^(E - FLINT_BITS dn) (normalized mantissa, u > 0,
   u != 1) to a relative accuracy of 2^-(func_prec + extra), into fix */
void _nfloat_log_fix(_nfloat_fix_t fix, nn_srcptr d, slong dn, slong E, slong extra, gr_ctx_t ctx);


/* For (y, len) > 0 known within err units of 2^ue, sets (t, wn) to its
   leading limbs normalized (top bit set), so that Y = T 2^(*te) within the
   returned number of units of 2^(*te) (shifted-out bits count as one more
   unit). Returns UWORD_MAX if the scaled bound would exceed 2^(FLINT_BITS - 3). */
ulong _nfloat_elem_normalize(nn_ptr t, slong * te, nn_srcptr y, slong len, slong ue, ulong err, slong wn);

/* Sets (q, w + 1) to floor(A B^w / D) or one more for normalized w-limb
   integers A and D (top bits set) known within ea resp. ed units; the
   quotient lies in (B^w / 2, 2 B^w) and the returned bound is in units
   of q. */
ulong _nfloat_elem_div(nn_ptr q, nn_srcptr a, ulong ea, nn_srcptr d, ulong ed, slong w);

/* A positive approximate real number (d, w) 2^(exp - FLINT_BITS w) with a
   normalized w-limb mantissa (top bit set) known within err units of its
   last limb (err == UWORD_MAX means that the bound overflowed). Used to
   compose a few operations with rigorous relative error bounds. */
typedef struct
{
    ulong d[NFLOAT_ELEM_MAX_LIMBS + 2];
    slong exp;
    ulong err;
}
_nfloat_rf_struct;

typedef _nfloat_rf_struct _nfloat_rf_t[1];

/* r = (y, len) 2^ue within err units */
void _nfloat_rf_set_mpn(_nfloat_rf_t r, nn_srcptr y, slong len, slong ue, ulong err, slong w);
/* r = a * b, a / b, a + b, |a - b|; the last returns the sign of a - b
   (0 for positive) and fails with err = UWORD_MAX for total cancellation */
void _nfloat_rf_mul(_nfloat_rf_t r, const _nfloat_rf_t a, const _nfloat_rf_t b, slong w);
void _nfloat_rf_div(_nfloat_rf_t r, const _nfloat_rf_t a, const _nfloat_rf_t b, slong w);
void _nfloat_rf_inv(_nfloat_rf_t r, const _nfloat_rf_t a, slong w);
void _nfloat_rf_add(_nfloat_rf_t r, const _nfloat_rf_t a, const _nfloat_rf_t b, slong w);
int _nfloat_rf_sub(_nfloat_rf_t r, const _nfloat_rf_t a, const _nfloat_rf_t b, slong w);
/* res = (-1)^sgnbit r; GR_UNABLE if the bound overflowed */
int _nfloat_set_rf(nfloat_ptr res, const _nfloat_rf_t r, int sgnbit, slong w, gr_ctx_t ctx);

/* Ball evaluation through mp_real, for arguments beyond the native
   ranges (and precisions beyond the temporaries).

   _nfloat_mp_real_view sets v = x (normal) exactly as an mp_real reading
   the limbs of buf (n + 1 limbs, which v must not outlive), for the
   mp_real functions taking const arguments; v must not be cleared or
   written.

   _nfloat_set_mp_real sets res = x 2^shift for a ball x: the midpoint in
   the default mode, the endpoint in the rounding direction otherwise.
   It returns 0 without writing res if the ball contains zero or its
   relative radius exceeds 2^-(func_prec + NFLOAT_ELEM_GUARD_BITS), else
   1 (with the status of the assignment in *status).

   _nfloat_eval_ball sets res = f(x) (or f(x, y)) by calling
   f(r, &shift, X, Y, op, prec) on the exact arguments at increasing
   precision prec, starting from at least min_prec bits, until the ball
   r 2^shift determines res.  The callback returns 1 when r 2^shift
   contains the value, 0 when prec did not suffice to evaluate it (a
   ball not certified to lie in the domain, or too wide for a division),
   or -1 when the function is undefined there (then res = nan).  Returns
   GR_UNABLE if the precision exceeds 16 prec + 100000 bits (or 2 min_prec
   bits, for min_prec up to NFLOAT_TRIG_BALL_MAX_EXP / 4 + 64). */
typedef int (* _nfloat_ball_func_t)(mp_real_t r, slong * shift, const mp_real_t x, const mp_real_t y, int op, slong prec);

/* the minimum precision of the trigonometric functions of x in the
   ball evaluation: mp_real_sin_cos_bits and mp_real_tan_bits reduce
   arguments up to 2^max(65536, 4 prec); arguments beyond
   2^NFLOAT_TRIG_BALL_MAX_EXP are not reduced (GR_UNABLE) */
#define NFLOAT_TRIG_BALL_MAX_EXP (WORD(1) << 22)

FLINT_FORCE_INLINE slong
_nfloat_trig_ball_min_prec(nfloat_srcptr x)
{
    slong e = NFLOAT_EXP(x);

    if (e > NFLOAT_TRIG_BALL_MAX_EXP)
        return WORD_MAX / 4;
    return (e > 65536) ? e / 4 + 64 : 0;
}

void _nfloat_mp_real_view(mp_real_t v, nn_ptr buf, nfloat_srcptr x, gr_ctx_t ctx);
int _nfloat_set_mp_real(nfloat_ptr res, const mp_real_t x, slong shift, int * status, gr_ctx_t ctx);
int _nfloat_eval_ball(nfloat_ptr res, _nfloat_ball_func_t f, nfloat_srcptr x, nfloat_srcptr y, int op, slong min_prec, gr_ctx_t ctx);

#endif
