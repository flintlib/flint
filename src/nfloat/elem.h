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

#include "nfloat.h"
#include "impl.h"
#include "../arf/elem.h"

/*
    Conventions: see arf/elem.h, which holds the shared definitions
    (fixed-point frames, the register fast paths, the limb counts and
    error bounds).  The functions here assemble the fixed-point results
    into nfloats: truncated in the default rounding mode, and in the
    floor and ceiling modes the endpoint of the error interval in the
    rounding direction, rounded outward, so that the result is a valid
    bound.
*/

/* Limbs of temporaries: target precision plus guard bits plus as many
   bits again for arguments close to zeros of the function. */
#define NFLOAT_ELEM_MAX_LIMBS (2 * NFLOAT_MAX_LIMBS + 8)

/* The largest left shift in normalizing a fixed-point value of small
   magnitude to a floating-point mantissa, which scales up its error count:
   that must stay below 2^(FLINT_BITS - 4) (the bound of the _nfloat_rf_t
   operations) for kernel errors below 2^14 ulps and a few operations, so
   40 bits on 64-bit machines and 8 on 32-bit machines; values with more
   leading zero bits are computed at more limbs. */
#define NFLOAT_ELEM_MAX_NORM_SHIFT (FLINT_BITS - 24)

/* the bits for the kernels (mp_real.h) for an absolute accuracy of
   2^-(p + extra), p the function precision of the context, and the
   fixed-point limbs for them */
FLINT_FORCE_INLINE slong
_nfloat_elem_bits(gr_ctx_t ctx, slong extra)
{
    return NFLOAT_CTX_FUNC_PREC(ctx) + _mp_real_elem_guard_bits(NFLOAT_CTX_FUNC_PREC(ctx)) + extra;
}

FLINT_FORCE_INLINE slong
_nfloat_elem_limbs(gr_ctx_t ctx, slong extra)
{
    return _mp_real_elem_limbs(_nfloat_elem_bits(ctx, extra));
}

/* Sets res to (-1)^sgnbit V where V is a real number with
   (Y - elo) 2^e <= V <= (Y + ehi) 2^e for the integer Y = (y, len) > 0
   (y must have room for len + 1 limbs; it is destroyed). */
int _nfloat_set_mpn_err(nfloat_ptr res, nn_ptr y, slong len, slong e, ulong elo, ulong ehi, int sgnbit, gr_ctx_t ctx);

/* the return value of a fast path that does not apply to its arguments
   (on 32-bit machines, of every fast path) */
#define NFLOAT_FAST_DECLINE (-1)

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

/* the reciprocal trigonometric functions and sinc by the fast paths
   (sin_cos.c) for one and two limbs, or NFLOAT_FAST_DECLINE */
#define NFLOAT_TRIG_COT 1
#define NFLOAT_TRIG_SEC 2
#define NFLOAT_TRIG_CSC 3
#define NFLOAT_TRIG_SINC 4
int _nfloat_trig_recip_fast(nfloat_ptr res, nfloat_srcptr x, int pi, int which, gr_ctx_t ctx);

/* inverse functions at one and two limbs (the compositions of
   arf/elem_arc_fast.c), for normal x in the domain (|x| != 1 where that
   is special) */
#define NFLOAT_ARC_ASIN ARF_ARC_ASIN
#define NFLOAT_ARC_ACOS ARF_ARC_ACOS
#define NFLOAT_ARC_ATAN ARF_ARC_ATAN
#define NFLOAT_ARC_ACOT ARF_ARC_ACOT
#define NFLOAT_ARC_ASEC ARF_ARC_ASEC
#define NFLOAT_ARC_ACSC ARF_ARC_ACSC
#define NFLOAT_ARC_ASINH ARF_ARC_ASINH
#define NFLOAT_ARC_ACOSH ARF_ARC_ACOSH
#define NFLOAT_ARC_ATANH ARF_ARC_ATANH
#define NFLOAT_ARC_ACOTH ARF_ARC_ACOTH
#define NFLOAT_ARC_ASECH ARF_ARC_ASECH
#define NFLOAT_ARC_ACSCH ARF_ARC_ACSCH
int _nfloat_arc_fast(nfloat_ptr res, nfloat_srcptr x, int f, int div_pi, gr_ctx_t ctx);

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
   resp. -2^e < delta <= 0 (dsgn = 1), e <= -2: +-1 in the default
   rounding mode (within the function precision when e is at most about
   minus it), and a valid bound otherwise (1 +- 2^e, or one ulp from 1
   for e <= -FLINT_BITS * nlimbs). */
int _nfloat_one_plus_tiny(nfloat_ptr res, int sgnbit, int dsgn, slong e, gr_ctx_t ctx);

/* Sets res to x (1 + delta) with -2^e <= delta <= 0 (side = -1) resp.
   0 <= delta <= 2^e (side = 1), e < -FLINT_BITS * nlimbs (a tiny
   relative perturbation): x in the default rounding mode, a valid bound
   otherwise. */
int _nfloat_set_tiny_rel_err(nfloat_ptr res, nfloat_srcptr x, int side, gr_ctx_t ctx);

/* res = (-1)^neg Y (recip = 0) or (-1)^neg / Y (recip = 1) for
   Y = (y, len) 2^ue within err units (the reciprocal by one division at
   the working precision; GR_UNABLE for Y = 0) */
int _nfloat_set_mpn_err_recip(nfloat_ptr res, nn_ptr y, slong len, slong ue, ulong err, int recip, int neg, gr_ctx_t ctx);

#define _nfloat_elem_normalize _mp_real_fp_normalize

/* Sets (q, w + 1) to floor(A B^w / D) or one more for normalized w-limb
   integers A and D (top bits set) known within ea resp. ed units; the
   quotient lies in (B^w / 2, 2 B^w) and the returned bound is in units
   of q. */
ulong _nfloat_elem_div(nn_ptr q, nn_srcptr a, ulong ea, nn_srcptr d, ulong ed, slong w);

/* Ball evaluation through mp_real, for arguments beyond the native
   ranges (and precisions beyond the temporaries).

   _nfloat_mp_real_view sets v = x (normal) exactly as an mp_real reading
   the limbs of buf (n + 1 limbs, which v must not outlive), for the
   mp_real functions taking const arguments; v must not be cleared or
   written.

   _nfloat_set_mp_real sets res = x 2^shift for a ball x: the midpoint in
   the default mode, the endpoint in the rounding direction otherwise.
   It returns 0 without writing res if the ball contains zero or its
   relative radius exceeds 2^-(func_prec + guard bits), else
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
   arguments up to 2^MP_REAL_TRIG_MAGLIM(prec); arguments beyond
   2^NFLOAT_TRIG_BALL_MAX_EXP are not reduced (GR_UNABLE) */
#define NFLOAT_TRIG_BALL_MAX_EXP (WORD(1) << 22)

FLINT_FORCE_INLINE slong
_nfloat_trig_ball_min_prec(nfloat_srcptr x)
{
    slong e = NFLOAT_EXP(x);

    if (e > NFLOAT_TRIG_BALL_MAX_EXP)
        return WORD_MAX / 4;
    return (e > MP_REAL_TRIG_MAGLIM(0)) ? e / 4 + 64 : 0;
}

void _nfloat_mp_real_view(mp_real_t v, nn_ptr buf, nfloat_srcptr x, gr_ctx_t ctx);
int _nfloat_set_mp_real(nfloat_ptr res, const mp_real_t x, slong shift, int * status, gr_ctx_t ctx);
int _nfloat_eval_ball(nfloat_ptr res, _nfloat_ball_func_t f, nfloat_srcptr x, nfloat_srcptr y, int op, slong min_prec, gr_ctx_t ctx);

#endif
