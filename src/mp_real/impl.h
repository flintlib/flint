/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Internal declarations of the mp_real module: the module's translation
   units include this as "impl.h", the test, tuning and profiling
   programs as "mp_real/impl.h". */

#ifndef MP_REAL_IMPL_H
#define MP_REAL_IMPL_H

#include <math.h>
#include <string.h>
#include "longlong.h"
#include "double_extras.h"
#include "mpn_extras.h"
#include "mp_real.h"

#ifdef __cplusplus
extern "C" {
#endif

/* fudge factor absorbing the rounding of the radius computations that
   go through doubles (division, roots) */
#define MP_REAL_EPS (1.0 + 1e-6)

/* ball internals ***********************************************************/
/* ERROR BOUNDS.  A ball's radius is an integer count err < B of ulps
   of its bottom limb.  During the composition of a bound the
   contributions (the radii of the operands scaled by magnitudes, the
   truncations) live at limb anchors that may differ, so they are
   carried as pairs

       e = (hi B + lo) B^a,   a 128-bit count at the limb anchor a,

   and combined rounding UP in limb arithmetic: a sum aligns the lower
   anchor to the higher one, dividing its count by B^d and rounding up
   (one unit when d >= 2), so that every term costs at most one unit
   of the higher anchor beyond its value; a product of a radius by a
   magnitude (a top limb plus one) is one umul_ppmm.  Installing a
   bound on a ball then converts it to units of the mantissa's bottom
   limb, truncating limbs when it reaches B units (the dropped material
   is worth one more unit) and padding zero limbs when it lies below
   one ulp.  The bounds of the divisions and roots are composed in
   doubles (mp_real_dbnd_t, anchored the same way, rounding up with the fudge
   factor MP_REAL_EPS) and converted at the end. */

/* powers of the limb radix as doubles (exact) */
#if FLINT_BITS == 64
#define MP_REAL_D_B 0x1p64
#define MP_REAL_D_BINV 0x1p-64
#define MP_REAL_D_SQRTB 0x1p32
#define MP_REAL_D_SQRTBINV 0x1p-32
#define MP_REAL_LGB 6
#else
#define MP_REAL_D_B 0x1p32
#define MP_REAL_D_BINV 0x1p-32
#define MP_REAL_D_SQRTB 0x1p16
#define MP_REAL_D_SQRTBINV 0x1p-16
#define MP_REAL_LGB 5
#endif
/* B^2, the cap of the two-limb counts (hi, lo) */
#define MP_REAL_D_B2 (MP_REAL_D_B * MP_REAL_D_B)

/* a 128-bit count at a limb anchor */
typedef struct { ulong hi, lo; slong a; } mp_real_bnd_t;

#define _mp_real_bnd_zero ((mp_real_bnd_t) { 0, 0, 0 })

FLINT_FORCE_INLINE int
_mp_real_bnd_is_zero(mp_real_bnd_t x)
{
    return (x.hi | x.lo) == 0;
}

/* the radius of x */
FLINT_FORCE_INLINE mp_real_bnd_t
_mp_real_bnd_of(const mp_real_t x)
{
    mp_real_bnd_t e;
    e.hi = 0;
    e.lo = x->err;
    e.a = (x->size == 0) ? x->exp : x->exp - x->size;
    return e;
}

/* (hi, lo) B^a */
FLINT_FORCE_INLINE mp_real_bnd_t
_mp_real_bnd(ulong hi, ulong lo, slong a)
{
    mp_real_bnd_t e;
    e.hi = hi;
    e.lo = lo;
    e.a = a;
    return e;
}

/* MAGNITUDES.  The magnitude of a ball with nonzero mantissa is
   bounded through the leading FLINT_BITS bits m of its mantissa (across
   a limb boundary), normalized with the top bit set by a shift s:

       m 2^-s B^(exp - 1) <= |x| < (m + 1) 2^-s B^(exp - 1),

   so the bounds overstate by at most 2^-63 relatively.  Bounding
   through the top limb alone, |x| < (top + 1) B^(exp - 1), overstates
   by a factor up to 2 for a top limb of 1 (every value in [1, 2),
   the working set of any iteration near 1) and by up to B in the
   other direction; with the radius held as a whole count of ulps
   that pessimism turns directly into truncated limbs. */
FLINT_FORCE_INLINE ulong
_mp_real_lead(const mp_real_t x, int * s)
{
    ulong m = x->d[x->size - 1];
    int c = flint_clz(m);
    if (c != 0)
    {
        m <<= c;
        if (x->size > 1)
            m |= x->d[x->size - 2] >> (FLINT_BITS - c);
    }
    *s = c;
    return m;
}

/* lower and upper bounds of |x| / B^(exp - 1) as doubles in [1, B] */
FLINT_FORCE_INLINE double
_mp_real_mag_lo(const mp_real_t x)
{
    int s;
    ulong m = _mp_real_lead(x, &s);
    return d_mul_2exp_inrange((double) m * (1.0 - 0x1p-52), -s);
}

FLINT_FORCE_INLINE double
_mp_real_mag_hi(const mp_real_t x)
{
    int s;
    ulong m = _mp_real_lead(x, &s);
    return d_mul_2exp_inrange((double) m * (1.0 + 0x1p-52) + 1.0, -s);
}

/* v |x| B^(a + 1 - exp) bounded: v (m + 1) 2^-s B^a, rounded up */
FLINT_FORCE_INLINE mp_real_bnd_t
_mp_real_bnd_mul_mag(ulong v, const mp_real_t x, slong a)
{
    mp_real_bnd_t e;
    int s;
    ulong m = _mp_real_lead(x, &s), r;

    umul_ppmm(e.hi, e.lo, v, m);
    add_ssaaaa(e.hi, e.lo, e.hi, e.lo, UWORD(0), v);
    if (s != 0)
    {
        r = e.lo & ((UWORD(1) << s) - 1);
        e.lo = (e.lo >> s) | (e.hi << (FLINT_BITS - s));
        e.hi >>= s;
        if (r != 0)
            add_ssaaaa(e.hi, e.lo, e.hi, e.lo, UWORD(0), UWORD(1));
    }
    e.a = a;
    return e;
}

/* x + y rounded up: the lower anchor is aligned to the higher one, its
   count divided by B^d rounding up (one unit for d >= 2); a carry out
   of 128 bits moves the anchor up a limb */
FLINT_FORCE_INLINE mp_real_bnd_t
_mp_real_bnd_add(mp_real_bnd_t x, mp_real_bnd_t y)
{
    mp_real_bnd_t r;
    slong d;
    ulong yhi, ylo, cy;

    if (_mp_real_bnd_is_zero(x)) return y;
    if (_mp_real_bnd_is_zero(y)) return x;

    if (x.a < y.a)
    {
        r = x; x = y; y = r;
    }
    d = x.a - y.a;
    if (d == 0)
    {
        yhi = y.hi;
        ylo = y.lo;
    }
    else if (d == 1)
    {
        ylo = y.hi + (y.lo != 0);
        yhi = (ylo == 0);            /* y.hi = B - 1 with a nonzero low */
    }
    else
    {
        ylo = 1;
        yhi = 0;
    }
    add_ssaaaa(r.hi, r.lo, x.hi, x.lo, yhi, ylo);
    cy = (r.hi < x.hi) || (r.hi == x.hi && r.lo < x.lo);
    r.a = x.a;
    if (cy)
    {
        /* the sum is B^2 + (r.hi, r.lo): rounded up by a limb */
        ylo = r.hi + (r.lo != 0);
        r.hi = 1 + (ylo < r.hi);
        r.lo = ylo;
        r.a++;
    }
    return r;
}

/* the count of e as a single limb, rounded up, with its anchor: e is
   below 2 B^2 by construction, so at most two limbs move */
FLINT_FORCE_INLINE ulong
_mp_real_bnd_reduce(mp_real_bnd_t * e)
{
    if (e->hi != 0)
    {
        ulong v = e->hi + (e->lo != 0);
        e->a++;
        if (v == 0)
        {
            v = 1;
            e->a++;
            /* if e->hi + 1 overflowed, e was in [B^2 - B + 1, B^2):
               one unit of B^(a + 2) covers it */
        }
        e->hi = 0;
        e->lo = v;
    }
    return e->lo;
}

/* --- the double-anchored bounds of the divisions and roots --- */

typedef struct { double v; slong a; } mp_real_dbnd_t;

/* v == 0 or 1 <= v < B^2, restored by rebasing the anchor (exact
   power-of-two scalings; a limb pushed below the retained part is
   rounded up to one unit) */
static inline void
_mp_real_dbnd_reduce(mp_real_dbnd_t * x)
{
    if (x->v >= MP_REAL_D_B2)
    {
        x->v = x->v * MP_REAL_D_BINV + 1.0;
        x->a++;
        while (x->v >= MP_REAL_D_B2)
        {
            x->v = x->v * MP_REAL_D_BINV + 1.0;
            x->a++;
        }
    }
    else if (x->v < 1.0 && x->v != 0.0)
    {
        do
        {
            x->v = x->v * MP_REAL_D_B;
            x->a--;
        } while (x->v < 1.0);
    }
}

static inline mp_real_dbnd_t
_mp_real_dbnd(double v, slong a)
{
    mp_real_dbnd_t e;
    e.v = v;
    e.a = a;
    _mp_real_dbnd_reduce(&e);
    return e;
}

/* x + y, rounded up: aligned to the higher anchor, the lower operand
   scaled by B^-d for d = 1, 2 (exact) and dropped for d >= 3, where
   it is below 2^-64 of the higher one and thus inside the fudge
   factor MP_REAL_EPS applied to the sum */
static inline mp_real_dbnd_t
_mp_real_dbnd_add(mp_real_dbnd_t x, mp_real_dbnd_t y)
{
    mp_real_dbnd_t r;
    slong d;

    if (x.v == 0.0) return y;
    if (y.v == 0.0) return x;

    if (x.a < y.a)
    {
        r = x; x = y; y = r;
    }

    d = x.a - y.a;
    if (d == 0)
        r.v = x.v + y.v;
    else if (d == 1)
        r.v = x.v + y.v * MP_REAL_D_BINV;
    else if (d == 2)
        r.v = x.v + y.v * (MP_REAL_D_BINV * MP_REAL_D_BINV);
    else
        r.v = x.v;
    r.v *= MP_REAL_EPS;
    r.a = x.a;
    if (r.v >= MP_REAL_D_B2)
    {
        r.v = r.v * MP_REAL_D_BINV + 1.0;
        r.a++;
    }
    return r;
}

/* a double-anchored bound as a limb count: v < B^2 is split at
   B, the low part rounded up */
static inline mp_real_bnd_t
_mp_real_bnd_of_fberr(mp_real_dbnd_t e)
{
    mp_real_bnd_t r;
    double hi, lo;

    _mp_real_dbnd_reduce(&e);
    if (e.v == 0.0)
        return _mp_real_bnd_zero;
    hi = floor(e.v * MP_REAL_D_BINV);
    lo = ceil(e.v - hi * MP_REAL_D_B);
    r.a = e.a;
    r.hi = (ulong) hi;
    if (lo >= MP_REAL_D_B)
    {
        r.lo = 0;
        r.hi++;
        if (r.hi == 0)
        {
            /* B^2: one unit two limbs up */
            r.lo = 1;
            r.a += 2;
        }
    }
    else
        r.lo = (ulong) lo;
    return r;
}

void _mp_real_norm_slow(mp_real_t x);
FLINT_FORCE_INLINE void
_mp_real_norm(mp_real_t x)
{
    /* the common case: nonzero top limb, and a nonzero bottom limb
       unless inexact (where low zero limbs are kept as padding) */
    if (x->size > 0 && x->d[x->size - 1] != 0 && (x->err != 0 || x->d[0] != 0))
        return;
    _mp_real_norm_slow(x);
}

void _mp_real_apply_bnd_slow(mp_real_t x, mp_real_bnd_t e);
FLINT_FORCE_INLINE void
_mp_real_apply_bnd(mp_real_t x, mp_real_bnd_t e)
{
    /* the common cases: an exact result, or a count below B already
       at the scale of the mantissa's bottom limb */
    if (e.hi == 0)
    {
        if (e.lo == 0)
        {
            x->err = 0;
            return;
        }
        if (x->size > 0 && e.a == x->exp - x->size)
        {
            x->err = e.lo;
            return;
        }
    }
    _mp_real_apply_bnd_slow(x, e);
}

static inline void
_mp_real_apply_err(mp_real_t x, mp_real_dbnd_t e)
{
    _mp_real_apply_bnd(x, _mp_real_bnd_of_fberr(e));
}

/* set x to the ball 0 +/- e */
static inline void
_mp_real_zero_bnd(mp_real_t x, mp_real_bnd_t e)
{
    ulong v = _mp_real_bnd_reduce(&e);
    x->size = 0;
    x->negative = 0;
    x->exp = e.a;
    x->err = v;
}

static inline void
_mp_real_zero_err(mp_real_t x, mp_real_dbnd_t e)
{
    _mp_real_zero_bnd(x, _mp_real_bnd_of_fberr(e));
}

/* ACCURATE LIMBS.  With errors normalized below 2^69 ulps, an
   inexact mantissa is accurate to all but its ~2 lowest limbs, so
   the mantissa length itself tracks accuracy; the target length of
   an output whose relative accuracy is limited by operand x is
   size + 1 when x is inexact (one guard limb; anything longer only
   computes noise) and unbounded when x is exact. */
/* THE NOISE FLOOR of an inexact ball is its bottom limb: the radius
   is a count of at least one ulp of that limb, so no ball holds limbs
   below it, and any limb an operation computes below the floor of its
   result is dropped again when the radius is installed.  (The
   double-anchored radius of an earlier design could sit far below one
   ulp, and operations kept a pad of six limbs under the floor; in the
   limb design the pad only cost a pass over the operand -- a sum whose
   other operand reached below an inexact one moved it down and back
   up -- and six limbs of every product.)

   A sum therefore takes nothing below an inexact operand's floor (the
   discarded tail of the other costs one unit there), and products and
   quotients of inexact operands keep size + 2 limbs: the truncation
   below the second guard limb is then 1/B of the inherited radius. */
static inline slong
_mp_real_noise_floor(const mp_real_t x)
{
    return x->exp - x->size;
}

#define MP_REAL_ACC_GUARD 2

static inline slong
_mp_real_acc(const mp_real_t x, slong n)
{
    if (x->err == 0)
        return WORD_MAX;
    return x->size + MP_REAL_ACC_GUARD;
}

/* rigorous OVERestimate of the relative radius err / |x| of a ball
   with nonzero mantissa: |x| >= top B^(exp - 1), so
   rel <= (err / top) B^(1 - size); the scaling window is clamped on
   the small side to a value exceeding anything the invariant
   err < B allows there, keeping the estimate one-sided */
static inline double
_mp_real_rel_bound(const mp_real_t x)
{
    slong k = 1 - x->size;

    if (x->err == 0)
        return 0.0;
    if (x->size == 0)
        return HUGE_VAL;
    if (FLINT_BITS * k > 900)
        return HUGE_VAL;
    if (FLINT_BITS * k < -900)
        return 0x1p-700;    /* true bound < 2^(FLINT_BITS - 900) */
    /* |x| >= top B^(exp - 1): through the top limb, as the operation
       bounds do (the blanket B^(exp - 1) overstates by up to a limb) */
    return d_mul_2exp_inrange((double) x->err / _mp_real_mag_lo(x),
        (int) (FLINT_BITS * k)) * MP_REAL_EPS;
}

/* v B / c rounded up as a count at the anchor below v's: v exact in
   a double below 2^53 (c within 2^-53, the quotient within 2^-52, so
   a factor 1 + 2^-50 and one unit up cover the rounding), else by an
   integer division of the rare wide radius */
FLINT_FORCE_INLINE mp_real_bnd_t
_mp_real_bnd_div_ui_B(ulong v, ulong c, slong a)
{
    mp_real_bnd_t e;
    if (v < (UWORD(1) << FLINT_MIN(53, FLINT_BITS - 1)))
    {
        double d = (double) v / (double) c * (1.0 + 0x1p-50), f;
        e.hi = (ulong) d;
        f = (d - (double) e.hi) * MP_REAL_D_B;
        e.lo = (f >= MP_REAL_D_B - 2.0) ? (ulong) (MP_REAL_D_B - 2.0) : (ulong) f;
        add_ssaaaa(e.hi, e.lo, e.hi, e.lo, UWORD(0), UWORD(2));
    }
    else
    {
        ulong q, r;
        q = v / c;
        r = v - q * c;
        udiv_qrnnd(e.lo, r, r, UWORD(0), c);
        e.hi = q;
        add_ssaaaa(e.hi, e.lo, e.hi, e.lo, UWORD(0), UWORD(1));
    }
    e.a = a;
    return e;
}

/* ***************************************************************************/

/* tuning cutoffs */
#ifndef MP_REAL_NEWTON_CUTOFF
#define MP_REAL_NEWTON_CUTOFF 512
#endif
#ifndef MP_REAL_SIN_COS_NOTAB_BURST_CUTOFF
#define MP_REAL_SIN_COS_NOTAB_BURST_CUTOFF 1200
#endif
#ifndef MP_REAL_SIN_COS_NOTAB_SQRT_NEWTON_CUTOFF
#define MP_REAL_SIN_COS_NOTAB_SQRT_NEWTON_CUTOFF 2000
#endif

/* from-scratch ball evaluations of the constants (cached and exported
   by const_cache.c; pi itself for the formulas that need it) */
void _mp_real_const_pi_compute(mp_real_t res, slong n);
void _mp_real_const_pi4_compute(mp_real_t res, slong n);
void _mp_real_const_log2_compute(mp_real_t res, slong n);
void _mp_real_const_euler_compute(mp_real_t res, slong n);
void _mp_real_const_e_compute(mp_real_t res, slong n);
void _mp_real_const_log10_compute(mp_real_t res, slong n);
void _mp_real_const_catalan_compute(mp_real_t res, slong n);
void _mp_real_const_zeta3_compute(mp_real_t res, slong n);
void _mp_real_const_zeta5_compute(mp_real_t res, slong n);
void _mp_real_const_gamma_1_3_compute(mp_real_t res, slong n);
void _mp_real_const_gamma_1_4_compute(mp_real_t res, slong n);
void _mp_real_const_2_div_pi_compute(mp_real_t res, slong n);

/* ball-valued workers of the logarithm and arctangent (newton.c,
   log_agm.c) for a fixed-point input */
void _mp_real_neglog_newton_ball(mp_real_t res, nn_srcptr x, slong n, int forward, slong N);
void _mp_real_atan_newton_ball(mp_real_t res, nn_srcptr x, slong n, int forward, slong N);
void _mp_real_neglog_agm_ball(mp_real_t res, nn_srcptr x, slong n, slong N);

/* per-size kernels, generated by dev/tune_mp_real.py (64-bit only),
   and the generic series behind the _rs functions */
#if FLINT_BITS == 64
void _mp_real_exp_opt_1(nn_ptr res, nn_srcptr x);
void _mp_real_exp_opt_2(nn_ptr res, nn_srcptr x);
void _mp_real_exp_opt_3(nn_ptr res, nn_srcptr x);
void _mp_real_exp_opt_4(nn_ptr res, nn_srcptr x);
void _mp_real_exp_opt_5(nn_ptr res, nn_srcptr x);
void _mp_real_exp_opt_6(nn_ptr res, nn_srcptr x);
void _mp_real_exp_opt_7(nn_ptr res, nn_srcptr x);
void _mp_real_atan_opt_1(nn_ptr res, nn_srcptr x);
void _mp_real_atan_opt_2(nn_ptr res, nn_srcptr x);
void _mp_real_atan_opt_3(nn_ptr res, nn_srcptr x);
void _mp_real_atan_opt_4(nn_ptr res, nn_srcptr x);
void _mp_real_atan_opt_5(nn_ptr res, nn_srcptr x);
void _mp_real_atan_opt_6(nn_ptr res, nn_srcptr x);
void _mp_real_atan_opt_7(nn_ptr res, nn_srcptr x);
void _mp_real_log1p_opt_1(nn_ptr res, nn_srcptr x);
void _mp_real_log1p_opt_2(nn_ptr res, nn_srcptr x);
void _mp_real_log1p_opt_3(nn_ptr res, nn_srcptr x);
void _mp_real_log1p_opt_4(nn_ptr res, nn_srcptr x);
void _mp_real_log1p_opt_5(nn_ptr res, nn_srcptr x);
void _mp_real_log1p_opt_6(nn_ptr res, nn_srcptr x);
void _mp_real_log1p_opt_7(nn_ptr res, nn_srcptr x);
void _mp_real_trig_opt_1(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan, nn_srcptr x);
void _mp_real_trig_opt_2(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan, nn_srcptr x);
void _mp_real_trig_opt_3(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan, nn_srcptr x);
void _mp_real_trig_opt_4(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan, nn_srcptr x);
void _mp_real_trig_opt_5(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan, nn_srcptr x);
void _mp_real_trig_opt_6(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan, nn_srcptr x);
void _mp_real_trig_opt_7(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan, nn_srcptr x);
void _mp_real_trig_opt_8(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan, nn_srcptr x);
void _mp_real_trig_opt_9(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan, nn_srcptr x);
void _mp_real_trig_opt_10(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan, nn_srcptr x);
void _mp_real_trig_opt_11(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan, nn_srcptr x);
void _mp_real_trig_opt_12(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan, nn_srcptr x);
void _mp_real_sin_cos_opt_1(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void _mp_real_sin_cos_opt_2(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void _mp_real_sin_cos_opt_3(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void _mp_real_sin_cos_opt_4(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void _mp_real_sin_cos_opt_5(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void _mp_real_sin_cos_opt_6(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void _mp_real_sin_cos_opt_7(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void _mp_real_sin_cos_opt_8(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void _mp_real_sin_cos_opt_9(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void _mp_real_sin_cos_opt_10(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void _mp_real_sin_cos_opt_11(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void _mp_real_sin_cos_opt_12(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void _mp_real_tan_opt_1(nn_ptr res, nn_srcptr x);
void _mp_real_tan_opt_2(nn_ptr res, nn_srcptr x);
void _mp_real_tan_opt_3(nn_ptr res, nn_srcptr x);
void _mp_real_tan_opt_4(nn_ptr res, nn_srcptr x);
void _mp_real_tan_opt_5(nn_ptr res, nn_srcptr x);
void _mp_real_tan_opt_6(nn_ptr res, nn_srcptr x);
void _mp_real_tan_opt_7(nn_ptr res, nn_srcptr x);
void _mp_real_tan_opt_8(nn_ptr res, nn_srcptr x);
void _mp_real_tan_opt_9(nn_ptr res, nn_srcptr x);
void _mp_real_tan_opt_10(nn_ptr res, nn_srcptr x);
void _mp_real_tan_opt_11(nn_ptr res, nn_srcptr x);
void _mp_real_tan_opt_12(nn_ptr res, nn_srcptr x);
#endif

/* shared parts of the bitwise functions and the kernels */
void _mp_real_exp_recon(nn_ptr y, nn_ptr sh, slong ylen, const slong * used, slong j, slong num);
void _mp_real_tan_halfangle_mid(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_srcptr x, slong n, int r, void (*series)(nn_ptr, nn_srcptr));
/* tan, atan, atanh, sin or 1 - cos of t < 2^-r (r >= 8) at n limbs,
   (res, n), by the tapered Horner scheme with static coefficient
   tables (series_tapered.c), error below 6 ulps; the cost falls with
   r */
#define MP_REAL_SERIES_TAN 0
#define MP_REAL_SERIES_ATAN 1
#define MP_REAL_SERIES_ATANH 2
#define MP_REAL_SERIES_SIN 3
#define MP_REAL_SERIES_COS 4
#define MP_REAL_SERIES_SINH 5
#define MP_REAL_SERIES_COSH 6
#define MP_REAL_SERIES_EXP 7
#define MP_REAL_SERIES_EXP_NEG 8
/* tapered rectangular splitting (series_rs.c): exp, exp(-x) (res, n + 1), sin,
   sinh, atan, atanh, 1 - cos, cosh - 1 (res, n) at x < 2^-r, r >= 8;
   res must not alias x */
void _mp_real_series_rs(nn_ptr res, nn_srcptr x, slong n, flint_bitcnt_t r, int func);
void _mp_real_series_rs_sin_cos(nn_ptr ysin, nn_ptr yg, nn_srcptr x, slong n, flint_bitcnt_t r, int hyperbolic);
void _mp_real_series_tapered(nn_ptr res, nn_srcptr x, slong n, flint_bitcnt_t r, int func);
slong _mp_real_series_tapered_levels(slong n, flint_bitcnt_t r, int func);
slong _mp_real_series_tapered_horner(nn_ptr V, nn_srcptr z, slong n, flint_bitcnt_t r, int func, slong j0, slong N, slong g);
/* tan by the common-denominator chunks of series_rs.c and the tapered
   Horner tail, for r >= MP_REAL_SERIES_TAN_RMIN and n up to
   MP_REAL_SERIES_TAN_NMAX, beyond as long as the chunks cover every term
   (_mp_real_series_rs_tan_ok); within 3 ulps */
void _mp_real_series_rs_tan(nn_ptr res, nn_srcptr x, slong n, flint_bitcnt_t r);
int _mp_real_series_rs_tan_ok(slong n, flint_bitcnt_t r);

/* the static tables (elem_tables.c) cover n <= nmax, r >= rmin */
#include "elem_tables.h"

FLINT_FORCE_INLINE slong
_mp_real_series_tapered_nmax(int func)
{
    return (func == MP_REAL_SERIES_TAN) ? MP_REAL_SERIES_TAN_NMAX
         : (func <= MP_REAL_SERIES_ATANH) ? MP_REAL_SERIES_ATAN_NMAX
         : (func == MP_REAL_SERIES_SIN) ? MP_REAL_SERIES_SIN_NMAX
         : MP_REAL_SERIES_COS_NMAX;
}

FLINT_FORCE_INLINE slong
_mp_real_series_tapered_rmin(int func)
{
    return (func == MP_REAL_SERIES_TAN) ? MP_REAL_SERIES_TAN_RMIN
         : (func <= MP_REAL_SERIES_ATANH) ? MP_REAL_SERIES_ATAN_RMIN
         : (func == MP_REAL_SERIES_SIN) ? MP_REAL_SERIES_SIN_RMIN
         : MP_REAL_SERIES_COS_RMIN;
}

/* The reduced series of atan, atanh (func) resp. sin and 1 - cos at
   x < 2^-r, r >= 8, within 6 ulps: tapered Horner (about N/3 full
   products for N terms) below MP_REAL_SERIES_RS_MIN_* limbs, where it
   is faster, and tapered rectangular splitting (about 2 sqrt(N/3) full
   products plus scalar work) from there on (measured crossovers: 10-12
   limbs for atan/atanh at r = 10-16, 14 at r = 24-32, 16-18 for sin and
   cos at r = 24-64). */
#define MP_REAL_SERIES_RS_MIN_ATAN(r) (((r) >= 20) ? 14 : 12)
#define MP_REAL_SERIES_RS_MIN_SIN_COS 17

FLINT_FORCE_INLINE void
_mp_real_series_atan(nn_ptr res, nn_srcptr x, slong n, flint_bitcnt_t r, int func)
{
    if (n < MP_REAL_SERIES_RS_MIN_ATAN(r) && (slong) r >= MP_REAL_SERIES_ATAN_RMIN)
        _mp_real_series_tapered(res, x, n, r, func);
    else
        _mp_real_series_rs(res, x, n, r, func);
}

FLINT_FORCE_INLINE void
_mp_real_series_sin_cos(nn_ptr ysin, nn_ptr yg, nn_srcptr x, slong n, flint_bitcnt_t r)
{
    if (n < MP_REAL_SERIES_RS_MIN_SIN_COS && (slong) r >= MP_REAL_SERIES_SIN_RMIN)
    {
        _mp_real_series_tapered(ysin, x, n, r, MP_REAL_SERIES_SIN);
        _mp_real_series_tapered(yg, x, n, r, MP_REAL_SERIES_COS);
    }
    else
        _mp_real_series_rs_sin_cos(ysin, yg, x, n, r, 0);
}

/* floor(c B^n) for c = pi/4, log 2, 2/pi, n limbs (without copying):
   the static tables up to MP_REAL_CONST_STATIC_N limbs, the per-thread
   cache beyond (valid until the next growth of that entry) */
nn_srcptr _mp_real_const_cached_ptr(int which, slong n);
#define MP_REAL_CONST_ID_PI4 0
#define MP_REAL_CONST_ID_LOG2 1
#define MP_REAL_CONST_ID_2_DIV_PI 10

FLINT_FORCE_INLINE nn_srcptr
_mp_real_const_ptr(int which, slong n)
{
    if (n <= MP_REAL_CONST_STATIC_N)
    {
        const ulong * t = (which == MP_REAL_CONST_ID_PI4) ? _mp_real_const_pi4_static
            : (which == MP_REAL_CONST_ID_LOG2) ? _mp_real_const_log2_static
            : _mp_real_const_2_div_pi_static;
        return t + MP_REAL_CONST_STATIC_N - n;
    }
    return _mp_real_const_cached_ptr(which, n);
}
/* the half-angle reconstruction's range for the tangent series with a
   Horner table (series_tapered.c, and the tail of series_rs.c) */
#if FLINT_BITS == 64
#define MP_REAL_TAN_TAPERED_MIN 13
#ifndef MP_REAL_TAN_TAPERED_MAX
#define MP_REAL_TAN_TAPERED_MAX 80
#endif
/* beyond, the chunked tangent series of series_rs.c alone (when its
   chunks cover every term, _mp_real_series_rs_tan_ok) */
#ifndef MP_REAL_TAN_RS_MAX
#define MP_REAL_TAN_RS_MAX 600
#endif
#else
#define MP_REAL_TAN_TAPERED_MIN 1
#define MP_REAL_TAN_TAPERED_MAX 0
#define MP_REAL_TAN_RS_MAX 0
#endif
/* sin, cos, tan by the tangent half-angle reconstruction, any output
   NULL; returns 0 if the size is not handled */
int _mp_real_tan_halfangle(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan, nn_srcptr x, slong n, int r);

/* the greedy table-subtraction reduction; used must have room for
   MP_REAL_BITWISE_REDUCE_USED_ALLOC(r) entries */
#define MP_REAL_BITWISE_REDUCE_USED_ALLOC(r) \
    ((r) + 2 * ((r) / FLINT_BITS) + 12)
slong _mp_real_bitwise_reduce(nn_ptr t, slong wn, int r, slong istart,
    nn_srcptr tab, slong tabn, slong * used);

/* Per-thread tables of L_i = log(1 + 2^-i) and A_i = atan(2^-i),
   i = 0..r, each entry nv value limbs plus a guard limb below.  The
   storage is thread-local and cannot be exported from a Windows DLL:
   code outside the library reads entries through _entry, while the
   module uses the storage directly or through the inline _tab
   accessors, which prefer the static prefixes for small tables. */
void _mp_real_exp_logs_ensure(slong nv, slong rc);
nn_srcptr _mp_real_exp_logs_entry(slong i, slong n);
void _mp_real_exp_logs_clear(void);
slong _mp_real_exp_logs_max_index(void);
void _mp_real_atans_ensure(slong nv, slong rc);
nn_srcptr _mp_real_atans_entry(slong i, slong n);
void _mp_real_atans_clear(void);
slong _mp_real_atans_max_index(void);

extern FLINT_TLS_PREFIX nn_ptr _mp_real_exp_logs;
extern FLINT_TLS_PREFIX slong _mp_real_exp_logs_n;
extern FLINT_TLS_PREFIX slong _mp_real_exp_logs_r;
extern FLINT_TLS_PREFIX nn_ptr _mp_real_atans;
extern FLINT_TLS_PREFIX slong _mp_real_atans_n;
extern FLINT_TLS_PREFIX slong _mp_real_atans_r;

#define MP_REAL_STATIC_TAB_N 12
#define MP_REAL_STATIC_TAB_R 32
FLINT_DLL extern const ulong _mp_real_exp_logs_static[(MP_REAL_STATIC_TAB_R + 1) * MP_REAL_STATIC_TAB_N];
FLINT_DLL extern const ulong _mp_real_atans_static[(MP_REAL_STATIC_TAB_R + 1) * MP_REAL_STATIC_TAB_N];

FLINT_FORCE_INLINE nn_srcptr
_mp_real_exp_logs_tab(slong nv, slong rc, slong * nc)
{
    if (rc <= MP_REAL_STATIC_TAB_R && nv + 1 <= MP_REAL_STATIC_TAB_N)
    {
        *nc = MP_REAL_STATIC_TAB_N;
        return _mp_real_exp_logs_static;
    }
    _mp_real_exp_logs_ensure(nv, rc);
    *nc = _mp_real_exp_logs_n;
    return _mp_real_exp_logs;
}

FLINT_FORCE_INLINE nn_srcptr
_mp_real_atans_tab(slong nv, slong rc, slong * nc)
{
    if (rc <= MP_REAL_STATIC_TAB_R && nv + 1 <= MP_REAL_STATIC_TAB_N)
    {
        *nc = MP_REAL_STATIC_TAB_N;
        return _mp_real_atans_static;
    }
    _mp_real_atans_ensure(nv, rc);
    *nc = _mp_real_atans_n;
    return _mp_real_atans;
}

/* one-sided n-limb table values by binary splitting (floor or one
   below), i >= 1 */
void _mp_real_atan_2mexp_ui_bs(nn_ptr res, ulong i, slong n);
void _mp_real_log1p_2mexp_ui_bs(nn_ptr res, ulong i, slong n);

/* exact floors of table entries (tab_exact.c); entries within this
   many guard ulps of wrapping are recomputed exactly */
int _mp_real_tab_store_floor(nn_ptr e, const arb_t x, slong nc, slong prec);
void _mp_real_tab_entry_exact(nn_ptr e, int which, ulong i, slong nc);
#define MP_REAL_TAB_GUARD_SLACK UWORD(1024)

/* binary splitting of the reduced series (exp_sum_bs.c,
   sin_cos_sum_bs.c) */
slong _mp_real_exp_bs_num_terms(flint_bitcnt_t r, slong prec);
void _mp_real_exp_sum_bs_powtab(nn_ptr T, slong * tn, nn_ptr Q,
    slong * qn, slong * QE, nn_srcptr xp, slong xn, slong D, slong N);
void _mp_real_sin_cos_sum_bs_powtab(nn_ptr A, slong * an, slong * ae,
    nn_ptr B, slong * bn, slong * be, nn_ptr Q, slong * qn,
    slong * QE, nn_srcptr xp, slong xn, slong D, slong N, slong lmax);

/* the diophantine reductions: the relation-table descent and the
   per-thread caches of log p_j and 2 arg(pi_j) (nv fraction limbs and
   a units limb per entry, exact floors) */
typedef struct
{
    int gaussian;
    slong num;
    slong rows;
    const short * d;
    const double * epsilon;
}
mp_real_rel_static_struct;
FLINT_DLL extern const mp_real_rel_static_struct _mp_real_rel_static[];
FLINT_DLL extern const slong _mp_real_rel_static_num;

slong _mp_real_log_reduce(slong * rel, const mp_real_rel_struct * tab,
    nn_srcptr x, slong wr, double max_weight, double eps_min,
    nn_srcptr alpha, slong stride);
void _mp_real_log_dot(nn_ptr acc, nn_srcptr base, slong len,
    const slong * rel, slong num, nn_srcptr alpha, slong stride);
double _mp_real_signed_get_d(nn_srcptr a, slong len, nn_ptr tmp);

void _mp_real_log_primes_ensure(slong num, slong nv);
nn_srcptr _mp_real_log_primes_entry(slong j, slong nv);
void _mp_real_log_primes_clear(void);
slong _mp_real_log_primes_max_limbs(void);
void _mp_real_atan_gauss_ensure(slong num, slong nv);
nn_srcptr _mp_real_atan_gauss_entry(slong j, slong nv);
void _mp_real_atan_gauss_clear(void);

extern FLINT_TLS_PREFIX nn_ptr _mp_real_log_primes;
extern FLINT_TLS_PREFIX slong _mp_real_log_primes_n;
extern FLINT_TLS_PREFIX slong _mp_real_log_primes_num;
extern FLINT_TLS_PREFIX nn_ptr _mp_real_atan_gauss;
extern FLINT_TLS_PREFIX slong _mp_real_atan_gauss_n;
extern FLINT_TLS_PREFIX slong _mp_real_atan_gauss_num;

/* log p_j resp. 2 arg(pi_j), j < num, as balls, by the Machin-type
   sets (machin_bsplit.c); log(u/v) by Zuniga's series; exact floors of
   such values into nc-limb cache entries (0 if undetermined) */
void _mp_real_log_primes_vec(mp_real_struct * res, slong num, slong n);
void _mp_real_atan_gauss_vec(mp_real_struct * res, slong num, slong n);
void _mp_real_log_ratio_zuniga(mp_real_t res, const fmpz_t u, const fmpz_t v, slong n);
int _mp_real_store_floors(nn_ptr e, slong nc, mp_real_struct * v, slong num);

/* threads (parallel.c): run two functions, the second on a pool thread
   if one is free (returns 1 if so); run n tasks on up to n threads */
int _mp_real_parallel_pair(void (* f1)(void *), void * a1, void (* f2)(void *), void * a2);
void _mp_real_parallel_tasks(void (* f)(slong, void *), void * args, slong n);

/* T = T1 Q2 + P1 T2, Q = Q1 Q2, P = P1 P2 (need_p) in place, T2
   destroyed; par: on two threads.  A splitting forks its halves only
   while their exact values stay within MP_REAL_PAR_CAP times the
   working precision. */
void _mp_real_pqt_merge(mp_real_t P, mp_real_t Q, mp_real_t T, mp_real_t P2,
    mp_real_t Q2, mp_real_t T2, int need_p, slong n, int par);
#define MP_REAL_PAR_CAP 4.0

/* helpers of the elementary functions (sin_cos.c, exp.c, log.c,
   atan.c) ******************************************************************/

/* the radius of res += err B^anc, exactly */
FLINT_FORCE_INLINE void
_mp_real_elem_add_rad(mp_real_t res, ulong err, slong anc)
{
    _mp_real_apply_bnd(res, _mp_real_bnd_add(_mp_real_bnd_of(res),
        _mp_real_bnd(0, err, anc)));
}

/* the radius of res += |y| (midpoint magnitude plus radius), rounded
   up */
FLINT_FORCE_INLINE void
_mp_real_elem_add_mag(mp_real_t res, const mp_real_t y)
{
    mp_real_bnd_t e = _mp_real_bnd_of(y);
    if (y->size != 0)
        e = _mp_real_bnd_add(e, _mp_real_bnd_mul_mag(1, y, y->exp - 1));
    _mp_real_apply_bnd(res, _mp_real_bnd_add(_mp_real_bnd_of(res), e));
}

/* res->d holds (y, n + 1) with n fraction limbs: res = (-1)^neg y B^-n
   + [+- err ulps of B^-n], err >= 1 (low zero limbs are kept as the
   padding of an inexact value) */
FLINT_FORCE_INLINE void
_mp_real_elem_finish(mp_real_t res, slong n, ulong err, int neg)
{
    slong size = n + 1;

    while (size > 0 && res->d[size - 1] == 0)
        size--;

    if (size == 0)
    {
        mp_real_zero(res);
        res->exp = -n;
        _mp_real_elem_add_rad(res, err, -n);
        return;
    }

    res->size = size;
    res->exp = size - n;
    res->negative = neg;
    res->err = err;
}

/* |x| into the fixed-point frame (v, len) of n fraction limbs,
   truncated; returns 1 if limbs were dropped.  The caller guarantees
   that |x| < B^(len - n). */
FLINT_FORCE_INLINE int
_mp_real_elem_copy(nn_ptr v, slong len, slong n, const mp_real_t x)
{
    slong sh = x->exp - x->size + n;

    flint_mpn_zero(v, len);
    if (sh >= 0)
    {
        flint_mpn_copyi(v + sh, x->d, x->size);
        return 0;
    }
    if (x->size + sh > 0)
        flint_mpn_copyi(v, x->d - sh, x->size + sh);
    return 1;
}

/* the reduction parameters r of the specialized per-size kernels
   (exp_opt_<n>.c, log1p_opt_<n>.c, atan_opt_<n>.c, trig_opt_<n>.c on
   64-bit machines): the compile-time constants those files were
   emitted with by dev/tune_mp_real.py --pin; the kernels' default_r
   and the wrappers' guard-bit counts read them here */
#define MP_REAL_EXP_OPT_R_MAX 7
#define MP_REAL_EXP_OPT_R { 0, 12, 16, 16, 16, 16, 24, 32 }
#define MP_REAL_LOG1P_OPT_R_MAX 7
#define MP_REAL_LOG1P_OPT_R { 0, 16, 16, 10, 26, 31, 30, 25 }
#define MP_REAL_ATAN_OPT_R_MAX 7
#define MP_REAL_ATAN_OPT_R { 0, 4, 6, 22, 20, 18, 20, 19 }
#define MP_REAL_TRIG_OPT_R_MAX 12
#define MP_REAL_TRIG_OPT_R { 0, 4, 5, 9, 14, 15, 18, 16, 16, 16, 19, 23, 25 }

/* default_r without the call for the per-size range */
#if FLINT_BITS == 64
#define MP_REAL_DEFAULT_R_INLINE(name, NAME) \
FLINT_FORCE_INLINE int \
_mp_real_##name##_default_r_inline(slong n) \
{ \
    static const unsigned char tab[] = MP_REAL_##NAME##_OPT_R; \
    return (n <= MP_REAL_##NAME##_OPT_R_MAX) ? tab[n] \
        : _mp_real_##name##_bitwise_rs_default_r(n); \
}
#else
#define MP_REAL_DEFAULT_R_INLINE(name, NAME) \
FLINT_FORCE_INLINE int \
_mp_real_##name##_default_r_inline(slong n) \
{ \
    return _mp_real_##name##_bitwise_rs_default_r(n); \
}
#endif
MP_REAL_DEFAULT_R_INLINE(exp, EXP)
MP_REAL_DEFAULT_R_INLINE(log1p, LOG1P)
MP_REAL_DEFAULT_R_INLINE(atan, ATAN)
MP_REAL_DEFAULT_R_INLINE(trig, TRIG)

/* the per-size kernels called directly, as the bitwise functions'
   r = 0 dispatch does after its layers of checks: return 0 (nothing
   done) beyond the per-size range, else 1 with *err the bound those
   functions return (exp 9 r + 100, log1p 3 r + 64, atan 4 r + 64,
   sin/cos 6 r + 128) */
#if FLINT_BITS == 64
#define MP_REAL_OPT_CASE_1(name, k, res, x) \
    case k: _mp_real_##name##_opt_##k(res, x); break;
FLINT_FORCE_INLINE int
_mp_real_exp_opt(nn_ptr res, ulong * err, nn_srcptr x, slong n)
{
    switch (n)
    {
        MP_REAL_OPT_CASE_1(exp, 1, res, x) MP_REAL_OPT_CASE_1(exp, 2, res, x)
        MP_REAL_OPT_CASE_1(exp, 3, res, x) MP_REAL_OPT_CASE_1(exp, 4, res, x)
        MP_REAL_OPT_CASE_1(exp, 5, res, x) MP_REAL_OPT_CASE_1(exp, 6, res, x)
        MP_REAL_OPT_CASE_1(exp, 7, res, x)
        default: return 0;
    }
    *err = 9 * (ulong) _mp_real_exp_default_r_inline(n) + 100;
    return 1;
}

FLINT_FORCE_INLINE int
_mp_real_log1p_opt(nn_ptr res, ulong * err, nn_srcptr x, slong n)
{
    switch (n)
    {
        MP_REAL_OPT_CASE_1(log1p, 1, res, x) MP_REAL_OPT_CASE_1(log1p, 2, res, x)
        MP_REAL_OPT_CASE_1(log1p, 3, res, x) MP_REAL_OPT_CASE_1(log1p, 4, res, x)
        MP_REAL_OPT_CASE_1(log1p, 5, res, x) MP_REAL_OPT_CASE_1(log1p, 6, res, x)
        MP_REAL_OPT_CASE_1(log1p, 7, res, x)
        default: return 0;
    }
    *err = 3 * (ulong) _mp_real_log1p_default_r_inline(n) + 64;
    return 1;
}

FLINT_FORCE_INLINE int
_mp_real_atan_opt(nn_ptr res, ulong * err, nn_srcptr x, slong n)
{
    switch (n)
    {
        MP_REAL_OPT_CASE_1(atan, 1, res, x) MP_REAL_OPT_CASE_1(atan, 2, res, x)
        MP_REAL_OPT_CASE_1(atan, 3, res, x) MP_REAL_OPT_CASE_1(atan, 4, res, x)
        MP_REAL_OPT_CASE_1(atan, 5, res, x) MP_REAL_OPT_CASE_1(atan, 6, res, x)
        MP_REAL_OPT_CASE_1(atan, 7, res, x)
        default: return 0;
    }
    *err = 4 * (ulong) _mp_real_atan_default_r_inline(n) + 64;
    return 1;
}

#define MP_REAL_OPT_CASE_SC(k) \
    case k: _mp_real_trig_opt_##k(ysin, ycos, NULL, x); break;
FLINT_FORCE_INLINE int
_mp_real_sin_cos_opt(nn_ptr ysin, nn_ptr ycos, ulong * err, nn_srcptr x, slong n)
{
    switch (n)
    {
        MP_REAL_OPT_CASE_SC(1) MP_REAL_OPT_CASE_SC(2) MP_REAL_OPT_CASE_SC(3)
        MP_REAL_OPT_CASE_SC(4) MP_REAL_OPT_CASE_SC(5) MP_REAL_OPT_CASE_SC(6)
        MP_REAL_OPT_CASE_SC(7) MP_REAL_OPT_CASE_SC(8) MP_REAL_OPT_CASE_SC(9)
        MP_REAL_OPT_CASE_SC(10) MP_REAL_OPT_CASE_SC(11) MP_REAL_OPT_CASE_SC(12)
        default: return 0;
    }
    *err = 6 * (ulong) _mp_real_trig_default_r_inline(n) + 128;
    return 1;
}
#else
#define _mp_real_exp_opt(res, err, x, n) 0
#define _mp_real_log1p_opt(res, err, x, n) 0
#define _mp_real_atan_opt(res, err, x, n) 0
#define _mp_real_sin_cos_opt(ysin, ycos, err, x, n) 0
#endif

/* X[0, N) += q L[0, N) resp. X[0, N) -= q L[0, N), returning the carry
   resp. borrow limb, in longlong.h double-limb arithmetic (the sum
   L[i] q + cy + X[i] < B^2 never overflows two limbs); for a constant
   N the loop unrolls with everything in registers */
FLINT_FORCE_INLINE ulong
_mp_real_addmul_1_small(nn_ptr X, nn_srcptr L, slong N, ulong q)
{
    ulong cy = 0, hi, lo;
    slong i;

    for (i = 0; i < N; i++)
    {
        umul_ppmm(hi, lo, L[i], q);
        add_ssaaaa(hi, lo, hi, lo, 0, cy);
        add_ssaaaa(cy, X[i], hi, lo, 0, X[i]);
    }
    return cy;
}

FLINT_FORCE_INLINE ulong
_mp_real_submul_1_small(nn_ptr X, nn_srcptr L, slong N, ulong q)
{
    ulong cy = 0, hi, lo, t;
    slong i;

    for (i = 0; i < N; i++)
    {
        umul_ppmm(hi, lo, L[i], q);
        add_ssaaaa(hi, lo, hi, lo, 0, cy);
        sub_ddmmss(t, X[i], 0, X[i], hi, lo);
        cy = -t;
    }
    return cy;
}

/* res = X[0, N) +-= q L[0, N) with the register versions unrolled for
   N <= 6 and GMP's beyond */
#define MP_REAL_ADDMUL_1(res, X, L, N, q) \
    do { \
        switch (N) \
        { \
            case 1: (res) = _mp_real_addmul_1_small(X, L, 1, q); break; \
            case 2: (res) = _mp_real_addmul_1_small(X, L, 2, q); break; \
            case 3: (res) = _mp_real_addmul_1_small(X, L, 3, q); break; \
            case 4: (res) = _mp_real_addmul_1_small(X, L, 4, q); break; \
            case 5: (res) = _mp_real_addmul_1_small(X, L, 5, q); break; \
            case 6: (res) = _mp_real_addmul_1_small(X, L, 6, q); break; \
            default: (res) = mpn_addmul_1(X, L, N, q); break; \
        } \
    } while (0)

#define MP_REAL_SUBMUL_1(res, X, L, N, q) \
    do { \
        switch (N) \
        { \
            case 1: (res) = _mp_real_submul_1_small(X, L, 1, q); break; \
            case 2: (res) = _mp_real_submul_1_small(X, L, 2, q); break; \
            case 3: (res) = _mp_real_submul_1_small(X, L, 3, q); break; \
            case 4: (res) = _mp_real_submul_1_small(X, L, 4, q); break; \
            case 5: (res) = _mp_real_submul_1_small(X, L, 5, q); break; \
            case 6: (res) = _mp_real_submul_1_small(X, L, 6, q); break; \
            default: (res) = mpn_submul_1(X, L, N, q); break; \
        } \
    } while (0)

/* |x| < 1 as the fraction (v, n): read in place when x's top limb is
   the frame's top limb and x has at least n limbs (the common case),
   unless x shares its limbs with avoid1 or avoid2 (outputs the caller
   will write); else copied into buf (n limbs).  *trunc as for
   _mp_real_elem_copy. */
FLINT_FORCE_INLINE nn_srcptr
_mp_real_elem_frame(nn_ptr buf, slong n, const mp_real_t x,
    nn_srcptr avoid1, nn_srcptr avoid2, int * trunc)
{
    slong sh = x->exp - x->size + n;

    if (x->exp == 0 && sh <= 0 && x->d != avoid1 && x->d != avoid2)
    {
        *trunc = (sh < 0);
        return x->d - sh;
    }
    *trunc = _mp_real_elem_copy(buf, n, n, x);
    return buf;
}

/* leading zero bits of the fraction (v, n), WORD_MAX if zero */
FLINT_FORCE_INLINE slong
_mp_real_elem_lzb(nn_srcptr v, slong n)
{
    slong top = n - 1;
    while (top >= 0 && v[top] == 0)
        top--;
    return (top < 0) ? WORD_MAX : FLINT_BITS * (n - 1 - top) + flint_clz(v[top]);
}

/* (w, n) = s / (2 + s) resp. s / (2 - s) for the fraction (s, n)
   (0 < s < 1/2 for the minus form), rounded down or up by one ulp, by
   one approximate division (elem.c); w < 1.  With s within e ulps, w is
   within e + 1 ulps (the derivatives 2 / (2 +- s)^2 stay below 1 for
   s < 0.58). */
void _mp_real_elem_half_ratio(nn_ptr w, nn_srcptr s, slong n, int minus);

/* the radius of res += v B^a for a double v >= 0 (rounded up by the
   caller), through the double-anchored bounds */
FLINT_FORCE_INLINE void
_mp_real_elem_add_rad_d(mp_real_t res, double v, slong a)
{
    if (v != 0.0)
        _mp_real_apply_bnd(res, _mp_real_bnd_add(_mp_real_bnd_of(res),
            _mp_real_bnd_of_fberr(_mp_real_dbnd(v, a))));
}

/* an upper bound for rad(y) / B^(exp - 1) = err B^(1 - size): exact
   while FLINT_BITS (size - 1) <= 896 (then at least 2^-896), else the
   bound 2^-896 (err < B <= 2^64: the true value is below 2^(64 - 960)
   resp. 2^(32 - 928) on 32-bit), so that no subnormal arises */
FLINT_FORCE_INLINE double
_mp_real_elem_rad_rel_top(const mp_real_t y)
{
    if (y->err == 0)
        return 0.0;
    if (FLINT_BITS * (y->size - 1) > 896)
        return 0x1p-896;
    return d_mul_2exp_inrange((double) y->err, (int) (-FLINT_BITS * (y->size - 1)));
}

/* v >= |y| + rad(y), as v B^(exp - 1) for y with a nonzero mantissa */
FLINT_FORCE_INLINE double
_mp_real_elem_mag_hi(const mp_real_t y)
{
    return (_mp_real_mag_hi(y) + _mp_real_elem_rad_rel_top(y)) * (1.0 + 0x1p-52);
}

/* v <= |y| - rad(y) (possibly <= 0), as v B^(exp - 1) */
FLINT_FORCE_INLINE double
_mp_real_elem_mag_lo(const mp_real_t y)
{
    return (_mp_real_mag_lo(y) - _mp_real_elem_rad_rel_top(y)) * (1.0 - 0x1p-52);
}

/* 2^-e for 0 <= e <= 1000, else 0.0: a threshold whose tiny values
   may as well be zero, without passing through subnormals */
FLINT_FORCE_INLINE double
_mp_real_d_2exp_neg_or_zero(slong e)
{
    FLINT_ASSERT(e >= 0);
    return (e <= 1000) ? d_mul_2exp_inrange(1.0, (int) -e) : 0.0;
}

/* v 2^(FLINT_BITS a) for 2^-64 <= v <= 2^64 when a is at most 1 (larger
   a is the caller's to exclude), rounded up: exact while
   FLINT_BITS a > -960 (then at least 2^(-64 - 928), normal, on either
   word size), else the bound 2^-896 >= 2^64 2^-960 >= v 2^(FLINT_BITS a),
   so that the result is always a normal double */
FLINT_FORCE_INLINE double
_mp_real_elem_scale_up(double v, slong a)
{
    FLINT_ASSERT(a <= 1);
    if (FLINT_BITS * a <= -960)
        return 0x1p-896;
    return d_mul_2exp_inrange(v, (int) (FLINT_BITS * a));
}

/* res = m (exact, nonzero) truncated to its top k limbs, with one ulp
   of radius if limbs were dropped; m may share res's limbs (the
   midpoint of an aliased input) */
FLINT_FORCE_INLINE void
_mp_real_elem_set_trunc(mp_real_t res, const mp_real_t m, slong k)
{
    slong size = m->size, e = m->exp;
    int neg = m->negative, dropped = 0;
    nn_srcptr src;
    slong j;

    k = FLINT_MIN(k, size);
    for (j = 0; j < size - k && !dropped; j++)
        dropped = (m->d[j] != 0);
    src = m->d + size - k;
    if (res->d == m->d)
        memmove(res->d, src, k * sizeof(ulong));
    else
    {
        mp_real_fit_length(res, k);
        flint_mpn_copyi(res->d, src, k);
    }
    res->size = k;
    res->exp = e;
    res->negative = neg;
    res->err = dropped;
}

/* the ball [0 +- 2^e] resp. [1 +- 2^e] */
FLINT_FORCE_INLINE void
_mp_real_elem_set_error(mp_real_t res, int one, slong e)
{
    if (one)
        mp_real_set_ui(res, 1);
    else
        mp_real_zero(res);
    mp_real_add_error_2exp_si(res, e);
}

#ifdef __cplusplus
}
#endif

#endif
