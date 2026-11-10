/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Internal definitions for the elementary functions of arb: the point
   evaluations (an arf argument into a ball) over the fixed-point kernels
   of mp_real and the register fast paths of arf/elem.h. */

#ifndef ARB_ELEM_H
#define ARB_ELEM_H

#include "arb.h"
#include "arf/elem.h"

/*
    Conventions: the kernels return a fixed-point frame (y, len) 2^e
    within an error count in units of 2^e (arf/elem.h, mp_real.h); the
    ball is the frame rounded to prec bits (ARB_RND) with the error and
    the rounding error as the radius (_arb_set_mpn_err).  The kernels run
    at ceil((prec + _mp_real_elem_guard_bits(prec)) / FLINT_BITS) limbs, as for
    nfloat, so that the radius stays within a few ulps.

    The ball functions of arb (arb_exp, arb_sin_cos, ...) bound their
    working precision by the accuracy of the input and handle wide balls
    by monotonicity or endpoint evaluation; the point evaluations here
    (arb_exp_arf, arb_log_arf, arb_atan_arf, _arb_sin_cos_arf, ...) do the
    midpoint.
*/

/* the bits for the kernels (mp_real.h) and their limbs, for prec bits */
FLINT_FORCE_INLINE slong
_arb_elem_bits(slong prec)
{
    return prec + _mp_real_elem_guard_bits(prec);
}

FLINT_FORCE_INLINE slong
_arb_elem_limbs(slong prec)
{
    return _mp_real_elem_limbs(_arb_elem_bits(prec));
}

/* res = (-1)^sgnbit V with (Y - elo) 2^e <= V <= (Y + ehi) 2^e for the
   integer Y = (y, len) (len + 1 limbs available): the midpoint Y 2^e
   rounded to nearest at prec bits, the radius max(elo, ehi) 2^e plus
   half an ulp (arb/elem.c).  The exponent of the unit is e, a word; a result
   with a bignum exponent is rounded with e = 0 and shifted by the
   caller (arb_mul_2exp_fmpz). */
void _arb_set_mpn_err(arb_t res, nn_ptr y, slong len, slong e, ulong elo, ulong ehi, int sgnbit, slong prec);

/* the short series of a small exact argument (arf/elem.c), as a ball:
   the caller adds the tail to the radius */
#define ARB_SMALL_SERIES_SIN ARF_SMALL_SERIES_SIN
#define ARB_SMALL_SERIES_TAN ARF_SMALL_SERIES_TAN
#define ARB_SMALL_SERIES_ATAN ARF_SMALL_SERIES_ATAN
#define ARB_SMALL_SERIES_SINH ARF_SMALL_SERIES_SINH
#define ARB_SMALL_SERIES_ASINH ARF_SMALL_SERIES_ASINH
#define ARB_SMALL_SERIES_ATANH ARF_SMALL_SERIES_ATANH
#define ARB_SMALL_SERIES_COS ARF_SMALL_SERIES_COS
#define ARB_SMALL_SERIES_COSH ARF_SMALL_SERIES_COSH
#define ARB_SMALL_SERIES_LOG1P_POS ARF_SMALL_SERIES_LOG1P_POS
#define ARB_SMALL_SERIES_LOG1P_NEG ARF_SMALL_SERIES_LOG1P_NEG
void _arb_small_arg_series(arb_t z, nn_srcptr d, slong xn, slong e, int sgnbit, int kind, slong prec);

/* the inverse functions asin, acos, asinh, acosh, atanh (and acot) by
   the register compositions of arf/elem_arc_fast.c at one and two
   limbs, with the radius propagated by a derivative bound (arb/elem.c);
   returns 0 when they do not apply (res untouched) */
int _arb_arc_fast(arb_t res, const arb_t x, int f, slong prec);

/* res = 2^k times the constant id (MP_REAL_CONST_ID_PI4, ... of
   mp_real/impl.h) from its tabulated floor at prec + 16 or more bits:
   the midpoint rounded to nearest (correctly: checked by
   arb/test/t-const_floor), the radius half an ulp (arb/elem.c); for prec
   below about FLINT_BITS (MP_REAL_CONST_STATIC_N - 1) bits */
void _arb_const_set_floor(arb_t res, int id, int k, slong prec);

/* res = the cached constant (a function of _ARB_DEF_CACHED_CONSTANT) at
   prec + 32 bits rounded to nearest at prec bits, its radius half an ulp
   plus the cached error (arb/elem.c) */
void _arb_const_round_cached(arb_t res, void (*cached)(arb_t, slong), slong prec);

#if FLINT_BITS == 64
/* the frames of the register fast paths (arf/elem.h): res = (-1)^sgnbit
   (y2, y1, y0) 2^(k - 2 FLINT_BITS) resp. (y3, y2, y1, y0) 2^(k - 3
   FLINT_BITS) within err units (64-bit only, as the paths) */
/* the radius err 2^eu plus, when inexact, the rounding error below
   2^(exp - prec - 1) of the midpoint (exp its exponent before the
   rounding, whose carry only coarsens the ulp): one mag_set, the sum
   of the two terms in units 2^eu (err < 2^20; the rounding term at
   most 2^40 units, else in units 2^(eu + d - 40) with a 2^-40 slack) */
FLINT_FORCE_INLINE void
_arb_small_set_rad(mag_t rad, ulong err, slong eu, int inexact, slong exp, slong prec)
{
    slong d;

    if (!inexact)
    {
        if (err == 0)
            mag_zero(rad);
        else
            mag_set_ui_2exp_si(rad, err, eu);
        return;
    }

    FLINT_ASSERT(err < (UWORD(1) << 20));
    d = exp - prec - 1 - eu;
    /* (d >= 0: an inexact midpoint has a bit of the frame, at or above
       the unit 2^eu, below its rounding point 2^(exp - prec)) */
    FLINT_ASSERT(d >= 0);
    if (d >= 40)
        mag_set_ui_2exp_si(rad, (UWORD(1) << 40) + 1 + ((d < 60) ? (err >> (d - 40)) : 0), eu + d - 40);
    else
        mag_set_ui_2exp_si(rad, err + (UWORD(1) << d), eu);
}

/* the three- resp. four-limb frames of the register paths: Y 2^(k - 2
   resp. 3 FLINT_BITS) within err units, the midpoint rounded to nearest
   (arf/elem.h), the radius as above */
FLINT_FORCE_INLINE void
_arb_small_set_1(arb_t res, ulong y2, ulong y1, ulong y0, slong k, ulong err, int sgnbit, slong prec)
{
    ulong hi, lo, rest;
    unsigned int c;
    slong exp;

    if (FLINT_UNLIKELY(prec > 2 * FLINT_BITS || err >= (UWORD(1) << 20) || (y2 | y1) == 0))
    {
        ulong t[4];
        t[0] = y0; t[1] = y1; t[2] = y2;
        _arb_set_mpn_err(res, t, 3, k - 2 * FLINT_BITS, err, err, sgnbit, prec);
        return;
    }

    if (y2 != 0)
    {
        c = flint_clz(y2);
        if (c == 0)
        {
            hi = y2; lo = y1; rest = y0;
        }
        else
        {
            hi = (y2 << c) | (y1 >> (FLINT_BITS - c));
            lo = (y1 << c) | (y0 >> (FLINT_BITS - c));
            rest = y0 << c;
        }
        exp = k + FLINT_BITS - c;
    }
    else
    {
        c = flint_clz(y1);
        if (c == 0)
        {
            hi = y1; lo = y0;
        }
        else
        {
            hi = (y1 << c) | (y0 >> (FLINT_BITS - c));
            lo = y0 << c;
        }
        rest = 0;
        exp = k - (slong) c;
    }

    {
        /* inexact when any bit below the rounding point is set */
        slong sb = 2 * FLINT_BITS - prec;
        int inexact = (rest != 0) || (sb > FLINT_BITS ? (lo != 0 || (hi << (2 * FLINT_BITS - sb)) != 0)
                                      : sb == FLINT_BITS ? lo != 0 : sb > 0 && (lo << (FLINT_BITS - sb)) != 0);
        _arf_set_frame2(arb_midref(res), hi, lo, rest >> (FLINT_BITS - 1), (rest << 1) != 0, exp, sgnbit, prec, ARF_RND_NEAR);
        _arb_small_set_rad(arb_radref(res), err, k - 2 * FLINT_BITS, inexact, exp, prec);
    }
}

FLINT_FORCE_INLINE void
_arb_small_set_2(arb_t res, ulong y3, ulong y2, ulong y1, ulong y0, slong k, ulong err, int sgnbit, slong prec)
{
    ulong hi, lo, r1, r0;
    unsigned int c;
    slong exp;

    if (FLINT_UNLIKELY(prec > 2 * FLINT_BITS || err >= (UWORD(1) << 20) || (y3 | y2) == 0))
    {
        ulong t[5];
        t[0] = y0; t[1] = y1; t[2] = y2; t[3] = y3;
        _arb_set_mpn_err(res, t, 4, k - 3 * FLINT_BITS, err, err, sgnbit, prec);
        return;
    }

    if (y3 != 0)
    {
        c = flint_clz(y3);
        if (c == 0)
        {
            hi = y3; lo = y2; r1 = y1; r0 = y0;
        }
        else
        {
            hi = (y3 << c) | (y2 >> (FLINT_BITS - c));
            lo = (y2 << c) | (y1 >> (FLINT_BITS - c));
            r1 = (y1 << c) | (y0 >> (FLINT_BITS - c));
            r0 = y0 << c;
        }
        exp = k + FLINT_BITS - c;
    }
    else
    {
        c = flint_clz(y2);
        if (c == 0)
        {
            hi = y2; lo = y1; r1 = y0;
        }
        else
        {
            hi = (y2 << c) | (y1 >> (FLINT_BITS - c));
            lo = (y1 << c) | (y0 >> (FLINT_BITS - c));
            r1 = y0 << c;
        }
        r0 = 0;
        exp = k - (slong) c;
    }

    {
        slong sb = 2 * FLINT_BITS - prec;
        int inexact = (r1 | r0) != 0 || (sb > FLINT_BITS ? (lo != 0 || (hi << (2 * FLINT_BITS - sb)) != 0)
                                         : sb == FLINT_BITS ? lo != 0 : sb > 0 && (lo << (FLINT_BITS - sb)) != 0);
        _arf_set_frame2(arb_midref(res), hi, lo, r1 >> (FLINT_BITS - 1), ((r1 << 1) | r0) != 0, exp, sgnbit, prec, ARF_RND_NEAR);
        _arb_small_set_rad(arb_radref(res), err, k - 3 * FLINT_BITS, inexact, exp, prec);
    }
}

#endif

/* log(1 + x) of the exact x, keeping the relative accuracy for small x
   (arb/log_arf.c) */
void _arb_log1p_arf(arb_t z, const arf_t x, slong prec);

/* sin and cos of x = [m +- r] (either output NULL) by the range on the
   interval, in double arithmetic (arb/sin_cos_wide.c): for wide balls */
void _arb_sin_cos_wide(arb_t s, arb_t c, const arf_t x, const mag_t r, slong prec);
void arb_sin_cos_wide(arb_t s, arb_t c, const arb_t x, slong prec);

/* the point evaluations: sin and cos (either NULL) or tan of the exact x
   (pi: of pi x), at prec bits (arb/sin_cos.c), with upper bounds for
   |sin x| and |cos x| into ms resp. mc (either NULL; the derivatives
   for the propagation of a radius, nearly free from the frames); |x|
   beyond 2^maglim gives [0 +- 1] (the tangent [0 +- inf]) */
void _arb_sin_cos_arf(arb_t rs, arb_t rc, arb_t rt, mag_t ms, mag_t mc, const arf_t x, int pi, slong prec, slong maglim);

/* sinh and cosh (either NULL), or tanh, of the exact x (arb/sinh_cosh.c),
   with an upper bound for cosh x into mc (or NULL); |x| beyond 2^maglim
   overflows ([0 +- inf], the hyperbolic tangent +-1 within a tiny error) */
void _arb_sinh_cosh_arf(arb_t rs, arb_t rc, arb_t rt, mag_t mc, const arf_t x, slong prec, slong maglim);

/* log |a + bi| (which = 0) or atan2(b, a) (which = 1) by the register
   paths shared with acf at one and two limbs (arb/elem.c); 0 when the
   path does not apply */
#if FLINT_BITS == 64
int _arb_log_hypot_atan2_small(arb_t res, const arb_t a, const arb_t b, int which, slong prec);
#endif

#endif
