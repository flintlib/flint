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

/* mp_real_hypgeom_series with the content removal (hypgeom_bsplit.c)
   forced on (gcd = 1) or off (gcd = 0), or chosen by the precision
   (gcd = -1, the default) */
void _mp_real_hypgeom_series(mp_real_t res,
    const mp_real_hypgeom_series_struct * s, slong n, int gcd);

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
void _mp_real_const_inv_log2_frac_compute(mp_real_t res, slong n);

/* ball-valued workers of the logarithm and arctangent (newton.c,
   log_agm.c) for a fixed-point input */
void _mp_real_neglog_newton_ball(mp_real_t res, nn_srcptr x, slong n, int forward, slong N);
void _mp_real_atan_newton_ball(mp_real_t res, nn_srcptr x, slong n, int forward, slong N);
void _mp_real_neglog_agm_ball(mp_real_t res, nn_srcptr x, slong n, slong N);

/* the dense-table per-size kernels at 4 to 10 limbs, generated by
   dev/tune_mp_real.py (64-bit only),
   and the generic series behind the _rs functions */
#if FLINT_BITS == 64
void _mp_real_exp_opt_4(nn_ptr res, nn_srcptr x);
void _mp_real_exp_opt_5(nn_ptr res, nn_srcptr x);
void _mp_real_exp_opt_6(nn_ptr res, nn_srcptr x);
void _mp_real_exp_opt_7(nn_ptr res, nn_srcptr x);
void _mp_real_exp_opt_8(nn_ptr res, nn_srcptr x);
void _mp_real_exp_opt_9(nn_ptr res, nn_srcptr x);
void _mp_real_exp_opt_10(nn_ptr res, nn_srcptr x);
void _mp_real_atan_opt_4(nn_ptr res, nn_srcptr x);
void _mp_real_atan_opt_5(nn_ptr res, nn_srcptr x);
void _mp_real_atan_opt_6(nn_ptr res, nn_srcptr x);
void _mp_real_atan_opt_7(nn_ptr res, nn_srcptr x);
void _mp_real_atan_opt_8(nn_ptr res, nn_srcptr x);
void _mp_real_atan_opt_9(nn_ptr res, nn_srcptr x);
void _mp_real_atan_opt_10(nn_ptr res, nn_srcptr x);
void _mp_real_log1p_opt_4(nn_ptr res, nn_srcptr x);
void _mp_real_log1p_opt_5(nn_ptr res, nn_srcptr x);
void _mp_real_log1p_opt_6(nn_ptr res, nn_srcptr x);
void _mp_real_log1p_opt_7(nn_ptr res, nn_srcptr x);
void _mp_real_log1p_opt_8(nn_ptr res, nn_srcptr x);
void _mp_real_log1p_opt_9(nn_ptr res, nn_srcptr x);
void _mp_real_log1p_opt_10(nn_ptr res, nn_srcptr x);
void _mp_real_sin_cos_opt_4(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void _mp_real_sin_cos_opt_5(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void _mp_real_sin_cos_opt_6(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void _mp_real_sin_cos_opt_7(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void _mp_real_sin_cos_opt_8(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void _mp_real_sin_cos_opt_9(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void _mp_real_sin_cos_opt_10(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
#endif

/* shared parts of the bitwise functions and the kernels */
void _mp_real_exp_recon(nn_ptr y, nn_ptr sh, slong ylen, const slong * used, slong j, slong num);
void _mp_real_tan_halfangle_mid(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_srcptr x, slong n, int r);
void _mp_real_tan_halfangle_finish(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_ptr T, nn_ptr D, slong n);
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
/* the series of the hardcoded tiers (series_hard.c, 64-bit): x < 2^-zb
   with zb = 64 (x's top limb zero) or 32.  exp: (res, xn + 1), N <= 21
   terms, xn = ceil(N zb / 64); atan (alt) and atanh: (res, xn), xn <= 35
   resp. 17; sin and cos (alt), sinh and cosh: (ysin, xn + 1),
   (ycos, xn + 1), either may be NULL, xn <= 22 resp. 10 */
void _mp_real_exp_series_hard(nn_ptr res, nn_srcptr x, slong N, int zb);
void _mp_real_atan_series_hard(nn_ptr res, nn_srcptr x, slong xn, int zb, int alt);
void _mp_real_sin_cos_series_hard(nn_ptr ysin, nn_ptr ycos, nn_srcptr x, slong xn, int zb, int alt);
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

/* floor(c B^n) for c = pi/4, log 2, 2/pi, 1/log 2 - 1, n limbs (without copying):
   the static tables up to MP_REAL_CONST_STATIC_N limbs, the shared
   cache beyond (the pointer valid until flint_cleanup_master) */
nn_srcptr _mp_real_const_cached_ptr(int which, slong n);
slong _mp_real_const_cached_limbs(int which);
#define MP_REAL_CONST_ID_PI4 0
#define MP_REAL_CONST_ID_LOG2 1
#define MP_REAL_CONST_ID_2_DIV_PI 10
#define MP_REAL_CONST_ID_INV_LOG2_FRAC 11

FLINT_FORCE_INLINE nn_srcptr
_mp_real_const_ptr(int which, slong n)
{
    if (n <= MP_REAL_CONST_STATIC_N)
    {
        const ulong * t = (which == MP_REAL_CONST_ID_PI4) ? _mp_real_const_pi4_static
            : (which == MP_REAL_CONST_ID_LOG2) ? _mp_real_const_log2_static
            : (which == MP_REAL_CONST_ID_2_DIV_PI) ? _mp_real_const_2_div_pi_static
            : _mp_real_const_inv_log2_frac_static;
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
/* (t, n) = D B^n / S within one ulp (flint_mpn_divapprox_fraction on
   D's significant limbs), for a fraction (D, n) below the divisor
   (S, n + 1), S[n] != 0: the residual quotients of the log1p and atan
   reductions */
FLINT_FORCE_INLINE void
_mp_real_div_fraction(nn_ptr t, nn_srcptr D, nn_srcptr S, slong n)
{
    slong k = n;

    while (k > 0 && D[k - 1] == 0)
        k--;
    flint_mpn_zero(t, n);
    if (k > 0)
        flint_mpn_divapprox_fraction(t, D, k, S, n + 1, n);
}
/* (y, n + 1) = tan((v, n)) by the sine and cosine kernels and a division,
   for n <= MP_REAL_DENSE_RT_MAX on 64-bit machines */
void _mp_real_tan_sc(nn_ptr y, ulong * err, nn_srcptr v, slong n);

/* the greedy table-subtraction reduction; used must have room for
   MP_REAL_BITWISE_REDUCE_USED_ALLOC(r) entries */
#define MP_REAL_BITWISE_REDUCE_USED_ALLOC(r) \
    ((r) + 2 * ((r) / FLINT_BITS) + 12)
slong _mp_real_bitwise_reduce(nn_ptr t, slong wn, int r, slong istart,
    nn_srcptr tab, slong tabn, slong * used);

/* bitwise reduction with signed digits, the full-precision arithmetic
   in radix 2^56 on vec4n lanes (bitwise_vec.c, 64-bit only), within the
   bounds of the corresponding _bitwise_rs functions, which they serve
   from the given number of limbs: (res, n + 1) = exp((x, n)),
   (res, n) = log(1 + (x, n)), (res, n) = atan((x, n)), and sin, cos, tan
   as by _mp_real_tan_halfangle, 32 <= r <= 64 n - 16 */
#if FLINT_BITS == 64
void _mp_real_exp_bitwise_vec(nn_ptr res, nn_srcptr x, slong n, int r);
void _mp_real_log1p_bitwise_vec(nn_ptr res, nn_srcptr x, slong n, int r);
void _mp_real_atan_bitwise_vec(nn_ptr res, nn_srcptr x, slong n, int r);
void _mp_real_tan_halfangle_vec(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_srcptr x, slong n, int r);
/* free the radix tables of the vector reductions (no concurrent users) */
void _mp_real_bitwise_vec_clear(void);
#ifndef MP_REAL_EXP_BITWISE_VEC_MIN
#define MP_REAL_EXP_BITWISE_VEC_MIN 4
#endif
/* below, the limb reduction's work on the growing significant length
   of its product keeps up with the vector reduction (both measured at
   their tuned r, within a few percent of each other from 12 limbs) */
#ifndef MP_REAL_LOG1P_BITWISE_VEC_MIN
#define MP_REAL_LOG1P_BITWISE_VEC_MIN 48
#endif
#ifndef MP_REAL_ATAN_BITWISE_VEC_MIN
#define MP_REAL_ATAN_BITWISE_VEC_MIN 4
#endif
#ifndef MP_REAL_TRIG_BITWISE_VEC_MIN
#define MP_REAL_TRIG_BITWISE_VEC_MIN 4
#endif
#endif

/* The shared tables of L_i = log(1 + 2^-i) and A_i = atan(2^-i),
   i = 0..r, each entry nv value limbs plus a guard limb below, and of
   the Machin-type constants log p_j and 2 arg(pi_j) (log_primes.c,
   atan_gauss.c): process-wide shared caches (thread_support.h), one
   read-only version at a time, built at the largest precision and
   index range requested so far -- the floors nest, so smaller requests
   read the top limbs of the entries -- and grown by replacement, the
   version's limb count by at least half (bounding the memory held by
   the retired versions) and the entries only as far as requested, the
   old ones copied.  A version lives until flint_cleanup_master(), so
   the returned pointers stay valid for the duration of a computation;
   the module reads the versions directly or through the inline _tab
   accessors, which prefer the static prefixes for small tables, while
   code outside the library uses the _entry functions. */
typedef struct
{
    slong nc;               /* limbs per entry */
    slong num;              /* entries (r + 1 for the bitwise tables) */
    nn_ptr tab;
}
_mp_real_tab_struct;

/* the shared cache protocol of the limb tables (shared_tab.c): a
   request, the sufficiency test and the version destructor for
   flint_shared_cache_ensure, and a new version grown from the current
   one (its entries copied when the limb count stays, the builder
   filling from *i_from) */
typedef struct
{
    slong nc;
    slong num;
}
_mp_real_tab_req_struct;

int _mp_real_tab_sufficient(const void * version, void * args);
void _mp_real_tab_free(void * version);
_mp_real_tab_struct * _mp_real_tab_new(const _mp_real_tab_struct * old,
    const _mp_real_tab_req_struct * R, slong * i_from);

/* acc (nc limbs) += or -= (if sub) t >> s, the bits below the bottom
   dropped (a subtraction then takes one more ulp), scratch of nc limbs */
void _mp_real_tab_addshift(nn_ptr acc, nn_srcptr t, slong nc, slong s,
    int sub, nn_ptr scratch);

/* the version covering nv value limbs and the indices 0..rc */
const _mp_real_tab_struct * _mp_real_exp_logs_ensure(slong nv, slong rc);
const _mp_real_tab_struct * _mp_real_atans_ensure(slong nv, slong rc);
/* the top n limbs of entry i of the current version */
nn_srcptr _mp_real_exp_logs_entry(slong i, slong n);
nn_srcptr _mp_real_atans_entry(slong i, slong n);
/* free the cache (no other thread may be using it: tuning programs) */
void _mp_real_exp_logs_clear(void);
void _mp_real_atans_clear(void);
/* the largest index the current version covers (0 if none) */
slong _mp_real_exp_logs_max_index(void);
slong _mp_real_atans_max_index(void);

/* geometric growth of a shared table's limb count: the request, or one
   and a half times the old version's when that is larger */
FLINT_FORCE_INLINE slong
_mp_real_tab_grow_nc(slong nc_old, slong nc)
{
    if (nc <= nc_old)
        return nc_old;
    return FLINT_MAX(nc, nc_old + nc_old / 2);
}

#define MP_REAL_STATIC_TAB_N 12
#define MP_REAL_STATIC_TAB_R 32
FLINT_DLL extern const ulong _mp_real_exp_logs_static[(MP_REAL_STATIC_TAB_R + 1) * MP_REAL_STATIC_TAB_N];
FLINT_DLL extern const ulong _mp_real_atans_static[(MP_REAL_STATIC_TAB_R + 1) * MP_REAL_STATIC_TAB_N];

FLINT_FORCE_INLINE nn_srcptr
_mp_real_exp_logs_tab(slong nv, slong rc, slong * nc)
{
    const _mp_real_tab_struct * T;

    if (rc <= MP_REAL_STATIC_TAB_R && nv + 1 <= MP_REAL_STATIC_TAB_N)
    {
        *nc = MP_REAL_STATIC_TAB_N;
        return _mp_real_exp_logs_static;
    }
    T = _mp_real_exp_logs_ensure(nv, rc);
    *nc = T->nc;
    return T->tab;
}

FLINT_FORCE_INLINE nn_srcptr
_mp_real_atans_tab(slong nv, slong rc, slong * nc)
{
    const _mp_real_tab_struct * T;

    if (rc <= MP_REAL_STATIC_TAB_R && nv + 1 <= MP_REAL_STATIC_TAB_N)
    {
        *nc = MP_REAL_STATIC_TAB_N;
        return _mp_real_atans_static;
    }
    T = _mp_real_atans_ensure(nv, rc);
    *nc = T->nc;
    return T->tab;
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

/* transformed power tables for these splittings (powtab_pre.c): for
   the table entries used by at least two merges, the transform of the
   entry computed once, for partners of up to ymax[i] limbs */
#ifndef MP_REAL_POWTAB_PRE_CUTOFF
#define MP_REAL_POWTAB_PRE_CUTOFF 150
#endif
/* entries whose products exceed this many limbs are not transformed:
   a transform takes about 2.8 doubles per limb of its longest product
   and a tree keeps about twice its largest, so this bounds the extra
   memory to some 50 MB per tree (at 10^7 digits the largest products
   are 400000 limbs, all kept) */
#ifndef MP_REAL_POWTAB_PRE_MAX_LIMBS
#define MP_REAL_POWTAB_PRE_MAX_LIMBS (WORD(1) << 20)
#endif

typedef struct
{
    void * pre;             /* the transforms, or NULL */
    double * storage;
    slong len;
}
mp_real_powtab_pre_struct;

/* for N terms split at fixed midpoints (ranges of 1 or 2 terms being
   leaves): the number of merges using each step xexp[i] and their
   longest right child */
void _mp_real_bs_step_uses(slong * uses, slong * rmax, const slong * xexp,
    slong len, slong N);
void _mp_real_powtab_pre_init(mp_real_powtab_pre_struct * T,
    nn_srcptr * xpow, const slong * xlen, const slong * uses,
    const slong * ymax, slong len);
void _mp_real_powtab_pre_clear(mp_real_powtab_pre_struct * T);
/* z = xpow[i] y (xlen[i] + yn limbs, z may overlap y) through the
   transform of entry i; returns 0, writing nothing, when the entry has
   none or the product is outside its range */
int _mp_real_powtab_pre_mul(nn_ptr z, const mp_real_powtab_pre_struct * T,
    slong i, nn_srcptr y, slong yn);

/* the diophantine reductions: the relation-table descent and the
   shared caches of log p_j and 2 arg(pi_j) (nv fraction limbs and
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

/* the shared tables of log p_j and 2 arg(pi_j), j < num (see the
   bitwise tables above): entries of nc limbs, a unit limb on top, nv
   value fraction limbs and a guard limb at the bottom.  The version
   covering num entries at nv value limbs; the top nv + 1 limbs of
   entry j of a version (the unit limb at index nv); the limbs the
   current version holds (0 if none); clearing (no concurrent users) */
const _mp_real_tab_struct * _mp_real_log_primes_ensure(slong num, slong nv);
const _mp_real_tab_struct * _mp_real_atan_gauss_ensure(slong num, slong nv);
void _mp_real_log_primes_clear(void);
void _mp_real_atan_gauss_clear(void);
slong _mp_real_log_primes_max_limbs(void);

FLINT_FORCE_INLINE nn_srcptr
_mp_real_machin_entry(const _mp_real_tab_struct * T, slong j, slong nv)
{
    FLINT_ASSERT(j >= 0 && j < T->num);
    FLINT_ASSERT(nv >= 1 && nv + 2 <= T->nc);
    return T->tab + j * T->nc + (T->nc - (nv + 1));
}

/* log p_j resp. 2 arg(pi_j), j < num, as balls, by the Machin-type
   sets (machin_bsplit.c); log(u/v) by Zuniga's series; exact floors of
   such values into nc-limb cache entries (0 if undetermined) */
void _mp_real_log_primes_vec(mp_real_struct * res, slong num, slong n);
void _mp_real_atan_gauss_vec(mp_real_struct * res, slong num, slong n);
void _mp_real_log_ratio_zuniga(mp_real_t res, const fmpz_t u, const fmpz_t v, slong n);
void _mp_real_atan_ratio_zuniga(mp_real_t res, const fmpz_t d, const fmpz_t t, slong n);
int _mp_real_store_floors(nn_ptr e, slong nc, mp_real_struct * v, slong num);

/* threads (parallel.c): run two functions, the second on a pool thread
   if one is free (returns 1 if so); run n tasks on up to n resp.
   max_workers threads; products of vectors (destroyed) by forking
   product trees */
int _mp_real_parallel_pair(void (* f1)(void *), void * a1, void (* f2)(void *), void * a2);
void _mp_real_parallel_tasks(void (* f)(slong, void *), void * args, slong n);
void _mp_real_parallel_tasks_max(void (* f)(slong, void *), void * args, slong n, slong max_workers);
void _mp_real_parallel_tasks_cost(void (* f)(slong, void *), void * args, slong n, slong max_workers, const double * cost);
void _mp_real_parallel_lanes(void (* f)(slong, slong, void *), void * args, slong n, slong nlanes, const double * cost);
slong _mp_real_parallel_lanes_count(slong n, slong max_workers);
void _mp_real_vec_prod(mp_real_t res, mp_real_struct * vec, slong len, slong n);
void _mp_real_vec_prod_complex(mp_real_t rre, mp_real_t rim, mp_real_struct * re, mp_real_struct * im, slong len, slong n);

/* T = T1 Q2 + P1 T2, Q = Q1 Q2, P = P1 P2 (need_p) in place, T2
   destroyed; par: on two threads.  A splitting forks its halves only
   while their exact values stay within MP_REAL_PAR_CAP times the
   working precision. */
void _mp_real_pqt_merge(mp_real_t P, mp_real_t Q, mp_real_t T, mp_real_t P2,
    mp_real_t Q2, mp_real_t T2, int need_p, slong n, int par);
#define MP_REAL_PAR_CAP 4.0

/* The bit-burst loops (exp_reduced.c, sin_cos_reduced.c) evaluate
   their slices as parallel tasks; the concurrency is capped to one (a
   serial loop, the parallelism then coming from the binary splitting
   inside each slice) from these many working limbs, the switches of
   arb_exp_arf_bb (10^9 bits) and arb_sin_cos_arf_bb (4 10^8 bits): a
   slice's splitting integers and power table come to several times
   the working precision, so running the slices side by side multiplies
   the peak memory by the thread count.  (A single switch at 2^22
   limbs, 2.7 10^8 bits, left both functions at parity with arb on 8
   threads at 3 10^8 bits.) */
#ifndef MP_REAL_EXP_BURST_SERIAL_LIMBS
#define MP_REAL_EXP_BURST_SERIAL_LIMBS (1000000000 / FLINT_BITS)
#endif
#ifndef MP_REAL_TRIG_BURST_SERIAL_LIMBS
#define MP_REAL_TRIG_BURST_SERIAL_LIMBS (400000000 / FLINT_BITS)
#endif
/* below this many working limbs the burst's tasks and the product
   trees stay on the calling thread (a wake-up costs microseconds) */
#ifndef MP_REAL_PAR_MIN_LIMBS
#define MP_REAL_PAR_MIN_LIMBS 64
#endif

/* verbose reports (verbose.c): _mp_real_log prints unconditionally
   (the caller checks the level), _mp_real_progress at most once per
   MP_REAL_VERBOSE_INTERVAL seconds across all threads; a phase logs
   its start at level 2 and its end at level 2 or when it took at
   least the interval at level 1 */
#define MP_REAL_VERBOSE_INTERVAL 2.0
typedef struct { double t0; const char * name; slong n; } mp_real_phase_t;
void _mp_real_log(const char * fmt, ...);
void _mp_real_progress(const char * fmt, ...);
double _mp_real_verbose_time(void);
void _mp_real_phase_start(mp_real_phase_t * ph, const char * name, slong n);
void _mp_real_phase_end(mp_real_phase_t * ph);
/* phases are reported from this many limbs on */
#define MP_REAL_VERBOSE_MIN_LIMBS 4096
#define MP_REAL_PHASE_START(ph, nm, n) \
    do { if (mp_real_get_verbose() && (n) >= MP_REAL_VERBOSE_MIN_LIMBS) \
        _mp_real_phase_start(&(ph), nm, n); else (ph).name = NULL; } while (0)
#define MP_REAL_PHASE_END(ph) \
    do { if ((ph).name != NULL) _mp_real_phase_end(&(ph)); } while (0)

/* helpers of the elementary functions (sin_cos.c, exp.c, log.c,
   atan.c) ******************************************************************/

/* (ys, n + 1), (yc, n + 1) = sin v, cos v for v in [0, 1) at n fraction
   limbs, both outputs required, by the kernel mp_real_sin_cos_bits
   uses (per-size, bitwise, diophantine or notab by n, or a reduced
   series for small v); err <= 6r + 128 with r <= 768 the bitwise
   parameter, at most 128 for the others (sin_cos.c) */
void _mp_real_sin_cos_kernel(nn_ptr ys, nn_ptr yc, ulong * err, nn_srcptr v, slong n);

/* the kernel limbs for a sine and cosine accurate to about 2^-p
   (absolute), including the guard bits for the kernel's bound and a
   reduction error of a few ulps (sin_cos.c) */
slong _mp_real_sin_cos_limbs(slong p);

/* sin x and cos x (rs, rc, either may be NULL) from v in [0, 1) at n
   fraction limbs, for x = a pi/2 + s v (s = +-1) up to eps ulps of
   B^-n, with the sine negated for xneg; v must not alias the outputs'
   limbs (sin_cos.c) */
void _mp_real_sin_cos_eval(mp_real_t rs, mp_real_t rc, nn_srcptr v, slong n,
    ulong eps, int a, int s, int xneg);

/* reduction of an exact x with |x| >= 1 (x->exp >= 1) mod pi/2:
   |x| = a pi/2 + s v with v in [0, pi/4 (1 + 2/B)], (W, n) = v within
   3 ulps; returns a + 4 [s < 0] (sin_cos.c) */
int _mp_real_trig_reduce(nn_ptr W, const mp_real_t x, slong n);

/* the exact reduction of an exact x mod 2 for the functions of pi x:
   |x| = a/2 + s t, t in [0, 1/4], returning a + 4 [s < 0], with
   t = (T, *tl) B^(*te - *tl), T's top limb nonzero (*tl = 0 for t = 0),
   *te <= 0; T has room for x->size limbs.  And (V, N) = (pi/4) 4t B^N
   (= pi t) at N fraction limbs, below the exact value by less than 3
   ulps (trig_pi.c) */
int _mp_real_trig_pi_reduce(nn_ptr T, slong * tl, slong * te, const mp_real_t x);
void _mp_real_trig_pi_v_fixed(nn_ptr V, nn_srcptr T, slong tl, slong te, slong N);

/* (y, n + 1) = tan v for v in [0, 1) at n fraction limbs (the tangent
   series or sin and 1 - cos of the reduced argument and a division for
   small v, else the per-size, bitwise or diophantine kernels, or sin,
   cos and a division beyond 65536 limbs); err <= 8r + 256 with r <= 768
   the bitwise parameter, 4 for the diophantine kernel, computed at run
   time (below 1000) for the divisions (tan.c) */
void _mp_real_tan_kernel(nn_ptr y, ulong * err, nn_srcptr v, slong n);

/* res = (-1)^neg / y for y = (Y, L) B^ey (L limbs, top zero limbs
   ignored), known within err B^ey and a further relative factor
   1 +- 2^rx (WORD_MIN: none), to about wn limbs by one Newton inverse
   (or, when y is known to a relative accuracy below 2^-12 only, by the
   reciprocal of its ball): returns 1, or 0 (res untouched) if the ball
   is not strictly positive.  Y must not alias res's limbs (tan.c) */
int _mp_real_inv_mpn(mp_real_t res, nn_srcptr Y, slong L, slong ey, ulong err,
    slong rx, slong wn, int neg);

/* res = 1/x at n limbs for a strictly positive ball x of any width
   (res may alias x): returns 0 (res untouched) if x is not strictly
   positive, else 1 (tan.c) */
int _mp_real_ball_inv_pos(mp_real_t res, const mp_real_t x, slong n);

/* (y, n + 1) with a units limb from a ball enclosing a value in [0, 2)
   (clamped at 0 from below); returns the bound in ulps of B^-n,
   UWORD_MAX when none fits a word (tan.c) */
ulong _mp_real_ball_get_fixed1(nn_ptr y, const mp_real_t x, slong n);

/* the kernel limbs for a tangent accurate to about 2^-p (absolute),
   including the guard bits for the kernel's bound and a reduction
   error of a few ulps (tan.c) */
slong _mp_real_tan_limbs(slong p);

/* the propagation of an argument radius R (a ball enclosing it) through
   the tangent, res enclosing tan m for the midpoint m: returns 0 if the
   ball [m +- R] may contain a pole, else adds a bound for
   |tan y - tan m| over the ball to the radius of res and returns 1
   (tan.c) */
int _mp_real_tan_add_arg_rad(mp_real_t res, const mp_real_t R);

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

/* the reduction parameters r of the per-size kernels that r = 0 selects
   in the _bitwise_rs functions on 64-bit machines (the error bounds of
   the form 9 r + 100 and so on, and the wrappers' guard-bit counts, read
   them here): at n = 1 .. 3 the table-driven kernels of small.h, whose
   remainders lie below 2^-r; at n = 4 .. 10 the dense-table kernels
   (exp_opt_<n>.c, log1p_opt_<n>.c, atan_opt_<n>.c, sin_cos_opt_<n>.c),
   the compile-time constants dev/tune_mp_real.py --pin emitted them
   with, r = 7 L - 1 the reduced radius after L levels of seven bits.
   All of them stay far below their bounds of this form.  Two values of
   r select the hardcoded series of exp_rs.c and trig_rs.c instead of a
   hand series (see dev/gen_mp_real_dense_bodies.py): r = 32 after five
   levels, for the tier of 2^-32, and r = 63 after nine levels and a top
   level on the rest of the top limb, for the tier of 2^-64.  tan takes
   the sine and cosine kernels and a division (_mp_real_tan_sc). */
#define MP_REAL_EXP_OPT_R_MAX 10
#define MP_REAL_EXP_OPT_R { 0, 14, 14, 21, 34, 41, 48, 55, 55, 63, 63 }
#define MP_REAL_LOG1P_OPT_R_MAX 10
#define MP_REAL_LOG1P_OPT_R { 0, 20, 20, 27, 34, 48, 48, 55, 55, 62, 62 }
#define MP_REAL_ATAN_OPT_R_MAX 10
#define MP_REAL_ATAN_OPT_R { 0, 9, 9, 9, 20, 27, 34, 41, 41, 41, 41 }
#define MP_REAL_TRIG_OPT_R_MAX 10
#define MP_REAL_TRIG_OPT_R { 0, 14, 14, 14, 34, 34, 41, 41, 32, 32, 32 }
/* the dense-table sin and cos kernels (sin_cos_opt_<n>.c) */
#define MP_REAL_SIN_COS_OPT_MIN 4
#define MP_REAL_SIN_COS_OPT_MAX 10
#define MP_REAL_SIN_COS_OPT_R { 0, 0, 0, 0, 34, 34, 41, 41, 32, 32, 32 }

/* the dense-table kernels beyond the static tables (dense_rt.c, 64-bit
   only): MP_REAL_DENSE_W + 1 .. MP_REAL_DENSE_RT_MAX fraction limbs,
   nine levels of seven bits and a last one on the rest of the top limb
   (the remainder below 2^-64, in rare cases below 2^-MP_REAL_DENSE_RT_R,
   for the hardcoded series of the 2^-64 tier) with the first two levels
   static (at MP_REAL_DENSE_RT_W limbs) and the others computed at run
   time as needed, shared by all threads; exp writes n + 1 limbs, log1p
   and atan n, sin and cos n + 1 each.  MP_REAL_DENSE_RT_MAX is at most
   MP_REAL_DENSE_RT_W */
#define MP_REAL_DENSE_RT_LEVELS 10
#define MP_REAL_DENSE_RT_R 63
#ifndef MP_REAL_DENSE_RT_MAX
#define MP_REAL_DENSE_RT_MAX 32
#endif
#if FLINT_BITS == 64
void _mp_real_exp_dense_rt(nn_ptr y, nn_srcptr v, slong n);
void _mp_real_log1p_dense_rt(nn_ptr y, nn_srcptr v, slong n);
void _mp_real_atan_dense_rt(nn_ptr y, nn_srcptr v, slong n);
void _mp_real_sin_cos_dense_rt(nn_ptr ysin, nn_ptr ycos, nn_srcptr v, slong n);
/* fill every entry of the run-time tables of a family that is not
   there yet (the kernels fill them as needed), in the tier serving the
   kernels at n limbs; returns its width (21 for n <= 21, else
   MP_REAL_DENSE_RT_MAX) */
#define MP_REAL_DENSE_RT_EXP 0
#define MP_REAL_DENSE_RT_LOG1P 1
#define MP_REAL_DENSE_RT_ATAN 2
slong _mp_real_dense_rt_prepare(int which, slong n);
slong _mp_real_dense_rt_allocated(void);
/* (tests) entry j of level k (3 <= k <= MP_REAL_DENSE_RT_LEVELS) of a
   family in the tier serving n limbs, filled if needed, into (r, w),
   and its argument x = (a[1] + a[0] / B) / B; returns w.  The entries
   are lower bounds within one ulp for log(1 + x), -log(1 - x) and
   atan(x).  The largest index at level k. */
slong _mp_real_dense_rt_entry(nn_ptr r, ulong * a, int which, slong n,
    slong k, ulong j);
slong _mp_real_dense_rt_jmax(int which, slong k);
#endif

/* default_r without the call for the per-size range */
#if FLINT_BITS == 64
#define MP_REAL_DEFAULT_R_INLINE(name, NAME) \
FLINT_FORCE_INLINE int \
_mp_real_##name##_default_r_inline(slong n) \
{ \
    static const unsigned char tab[] = MP_REAL_##NAME##_OPT_R; \
    return (n <= MP_REAL_##NAME##_OPT_R_MAX) ? tab[n] \
        : (MP_REAL_##NAME##_DENSE_RT && n <= MP_REAL_DENSE_RT_MAX) \
            ? MP_REAL_DENSE_RT_R \
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
/* which families the run-time dense kernels serve past the per-size
   kernels (sin and cos below, separately from tan) */
#define MP_REAL_EXP_DENSE_RT 1
#define MP_REAL_LOG1P_DENSE_RT 1
#define MP_REAL_ATAN_DENSE_RT 1
#define MP_REAL_TRIG_DENSE_RT 0
MP_REAL_DEFAULT_R_INLINE(exp, EXP)
MP_REAL_DEFAULT_R_INLINE(log1p, LOG1P)
MP_REAL_DEFAULT_R_INLINE(atan, ATAN)
MP_REAL_DEFAULT_R_INLINE(trig, TRIG)

/* the reduction parameter behind the sin/cos kernel's bound 6 r + 128
   (the dense sin_cos_opt_<n>.c on their range, else as for trig) */
FLINT_FORCE_INLINE int
_mp_real_sin_cos_default_r_inline(slong n)
{
#if FLINT_BITS == 64
    static const unsigned char tab[] = MP_REAL_SIN_COS_OPT_R;
    if (n >= MP_REAL_SIN_COS_OPT_MIN && n <= MP_REAL_SIN_COS_OPT_MAX)
        return tab[n];
    if (n > MP_REAL_SIN_COS_OPT_MAX && n <= MP_REAL_DENSE_RT_MAX)
        return MP_REAL_DENSE_RT_R;
#endif
    return _mp_real_trig_default_r_inline(n);
}

/* Short divisions for the two-limb kernels, computing the same quotients
   as flint_mpn_tdiv_qr with 3/2 divisions in registers (hardware 2/1
   divisions or a 3/2 inverse, according to FLINT_PREINVERT_LIMB_USE_NATIVE) */

#if FLINT_PREINVERT_LIMB_USE_NATIVE
#define _MP_REAL_UDIV_QR_3BY2(q, r1, r0, n2, n1, n0, d1, d0, dinv) \
    FLINT_MPN_UDIV_QR_3BY2_HW(q, r1, r0, n2, n1, n0, d1, d0)
#define _MP_REAL_PREINV1(d1, d0) UWORD(0)
#else
#define _MP_REAL_UDIV_QR_3BY2(q, r1, r0, n2, n1, n0, d1, d0, dinv) \
    FLINT_MPN_UDIV_QR_3BY2(q, r1, r0, n2, n1, n0, d1, d0, dinv)
#define _MP_REAL_PREINV1(d1, d0) flint_mpn_preinv1(d1, d0)
#endif

/* (q[1], q[0]) = floor((n3, n2, n1, n0) / (d1, d0)) for a normalized
   divisor (d1 >= B/2) with (n3, n2) < (d1, d0) */
FLINT_FORCE_INLINE void
_mp_real_divq_4_2(nn_ptr q, ulong n3, ulong n2, ulong n1, ulong n0, ulong d1, ulong d0)
{
    ulong q1, q0, r1, r0, s1, s0, dinv;

    dinv = _MP_REAL_PREINV1(d1, d0);
    (void) dinv;
    _MP_REAL_UDIV_QR_3BY2(q1, r1, r0, n3, n2, n1, d1, d0, dinv);
    _MP_REAL_UDIV_QR_3BY2(q0, s1, s0, r1, r0, n0, d1, d0, dinv);
    (void) s1;
    (void) s0;
    q[1] = q1;
    q[0] = q0;
}

/* (q[1], q[0]) = floor((y1, y0, 0, 0) / (s2, s1, s0)) for s2 != 0, given
   that the quotient is below B^2 (schoolbook with 3/2 divisions by the
   top two limbs of the normalized divisor, each quotient limb at most one
   too large, as in GMP's mpn_sbpi1_div_qr) */
FLINT_FORCE_INLINE void
_mp_real_divq_4_3z(nn_ptr q, ulong y1, ulong y0, ulong s2, ulong s1, ulong s0)
{
    ulong d2, d1, d0, u2, u1, u0, qq, r1, r0, p1, p0, dinv;
    unsigned int sh;
    int j;

    FLINT_ASSERT(s2 != 0);

    sh = flint_clz(s2);
    if (sh != 0)
    {
        d2 = (s2 << sh) | (s1 >> (FLINT_BITS - sh));
        d1 = (s1 << sh) | (s0 >> (FLINT_BITS - sh));
        d0 = s0 << sh;
        u2 = y1 >> (FLINT_BITS - sh);
        u1 = (y1 << sh) | (y0 >> (FLINT_BITS - sh));
        u0 = y0 << sh;
    }
    else
    {
        d2 = s2; d1 = s1; d0 = s0;
        u2 = 0; u1 = y1; u0 = y0;
    }

    /* the partial remainder (u2, u1, u0) < (d2, d1, d0), followed by zero
       limbs */
    dinv = _MP_REAL_PREINV1(d2, d1);
    (void) dinv;

    for (j = 1; j >= 0; j--)
    {
        if (FLINT_UNLIKELY(u2 == d2 && u1 == d1))
        {
            /* the quotient limb is B - 1 (u0 < d0): the remainder is
               (u2, u1, u0, 0) - (B - 1) D = D - (0, d0 - u0, 0) */
            /* (the difference goes through a temporary: sub_dddmmmsss
               may write u0 before it reads its other inputs) */
            ulong t = d0 - u0;
            qq = ~UWORD(0);
            sub_dddmmmsss(u2, u1, u0, d2, d1, d0, UWORD(0), t, UWORD(0));
        }
        else
        {
            int neg;

            _MP_REAL_UDIV_QR_3BY2(qq, r1, r0, u2, u1, u0, d2, d1, dinv);
            /* (r1, r0, 0) - qq d0, negative at most once */
            umul_ppmm(p1, p0, qq, d0);
            neg = (r1 == 0) & ((r0 < p1) | ((r0 == p1) & (p0 != 0)));
            sub_dddmmmsss(u2, u1, u0, r1, r0, UWORD(0), UWORD(0), p1, p0);
            if (FLINT_UNLIKELY(neg))
            {
                qq--;
                add_sssaaaaaa(u2, u1, u0, u2, u1, u0, d2, d1, d0);
            }
        }

        q[j] = qq;
    }
}

#include "small.h"
#include "dense_tables.h"
#include "dense.h"


/* the per-size kernels called directly, as the bitwise functions'
   r = 0 dispatch does after its layers of checks, and past them up to
   MP_REAL_DENSE_RT_MAX limbs the dense kernels with run-time tables
   (dense_rt.c): return 0 (nothing done) beyond that range, else 1 with
   *err the bound of the bitwise form (exp 9 r + 100, log1p 3 r + 64,
   atan 4 r + 64, sin/cos 6 r + 128, r = 7 L - 1 for the dense kernels),
   or at one to three limbs the table-driven kernels of small.h with
   their smaller bounds */
#if FLINT_BITS == 64
#define MP_REAL_OPT_CASE_1(name, k, res, x) \
    case k: _mp_real_##name##_opt_##k(res, x); break;
FLINT_FORCE_INLINE int
_mp_real_exp_opt(nn_ptr res, ulong * err, nn_srcptr x, slong n)
{
    /* the table-driven kernels at one to three limbs (small.h) */
    if (n == 1)
    {
        _mp_real_small_exp_1_mpn(res, err, x);
        return 1;
    }
    if (n == 2)
    {
        _mp_real_small_exp_2_mpn(res, err, x);
        return 1;
    }
    if (n == 3)
    {
        _mp_real_small_exp_3_mpn(res, err, x);
        return 1;
    }

    switch (n)
    {
        MP_REAL_OPT_CASE_1(exp, 4, res, x)
        MP_REAL_OPT_CASE_1(exp, 5, res, x) MP_REAL_OPT_CASE_1(exp, 6, res, x)
        MP_REAL_OPT_CASE_1(exp, 7, res, x)
        MP_REAL_OPT_CASE_1(exp, 8, res, x) MP_REAL_OPT_CASE_1(exp, 9, res, x)
        MP_REAL_OPT_CASE_1(exp, 10, res, x)
        default:
            if (n > MP_REAL_DENSE_W && n <= MP_REAL_DENSE_RT_MAX)
            {
                _mp_real_exp_dense_rt(res, x, n);
                break;
            }
            return 0;
    }
    *err = 9 * (ulong) _mp_real_exp_default_r_inline(n) + 100;
    return 1;
}

FLINT_FORCE_INLINE int
_mp_real_log1p_opt(nn_ptr res, ulong * err, nn_srcptr x, slong n)
{
    if (n == 1)
    {
        _mp_real_small_log1p_1_mpn(res, err, x);
        return 1;
    }
    if (n == 2)
    {
        _mp_real_small_log1p_2_mpn(res, err, x);
        return 1;
    }
    if (n == 3)
    {
        _mp_real_small_log1p_3_mpn(res, err, x);
        return 1;
    }

    switch (n)
    {
        MP_REAL_OPT_CASE_1(log1p, 4, res, x)
        MP_REAL_OPT_CASE_1(log1p, 5, res, x) MP_REAL_OPT_CASE_1(log1p, 6, res, x)
        MP_REAL_OPT_CASE_1(log1p, 7, res, x)
        MP_REAL_OPT_CASE_1(log1p, 8, res, x) MP_REAL_OPT_CASE_1(log1p, 9, res, x)
        MP_REAL_OPT_CASE_1(log1p, 10, res, x)
        default:
            if (n > MP_REAL_DENSE_W && n <= MP_REAL_DENSE_RT_MAX)
            {
                _mp_real_log1p_dense_rt(res, x, n);
                break;
            }
            return 0;
    }
    *err = 3 * (ulong) _mp_real_log1p_default_r_inline(n) + 64;
    return 1;
}

FLINT_FORCE_INLINE int
_mp_real_atan_opt(nn_ptr res, ulong * err, nn_srcptr x, slong n)
{
    if (n == 1)
    {
        _mp_real_small_atan_1_mpn(res, err, x);
        return 1;
    }
    if (n == 2)
    {
        _mp_real_small_atan_2_mpn(res, err, x);
        return 1;
    }
    if (n == 3)
    {
        _mp_real_small_atan_3_mpn(res, err, x);
        return 1;
    }

    switch (n)
    {
        MP_REAL_OPT_CASE_1(atan, 4, res, x)
        MP_REAL_OPT_CASE_1(atan, 5, res, x) MP_REAL_OPT_CASE_1(atan, 6, res, x)
        MP_REAL_OPT_CASE_1(atan, 7, res, x)
        MP_REAL_OPT_CASE_1(atan, 8, res, x) MP_REAL_OPT_CASE_1(atan, 9, res, x)
        MP_REAL_OPT_CASE_1(atan, 10, res, x)
        default:
            if (n > MP_REAL_DENSE_W && n <= MP_REAL_DENSE_RT_MAX)
            {
                _mp_real_atan_dense_rt(res, x, n);
                break;
            }
            return 0;
    }
    *err = 4 * (ulong) _mp_real_atan_default_r_inline(n) + 64;
    return 1;
}

FLINT_FORCE_INLINE int
_mp_real_sin_cos_opt(nn_ptr ysin, nn_ptr ycos, ulong * err, nn_srcptr x, slong n)
{
    if (n == 1)
    {
        _mp_real_small_sin_cos_1_mpn(ysin, ycos, err, x);
        return 1;
    }
    if (n == 2)
    {
        _mp_real_small_sin_cos_2_mpn(ysin, ycos, err, x);
        return 1;
    }
    if (n == 3)
    {
        _mp_real_small_sin_cos_3_mpn(ysin, ycos, err, x);
        return 1;
    }

    if (n >= MP_REAL_SIN_COS_OPT_MIN && n <= MP_REAL_SIN_COS_OPT_MAX)
    {
        switch (n)
        {
            case 4: _mp_real_sin_cos_opt_4(ysin, ycos, x); break;
            case 5: _mp_real_sin_cos_opt_5(ysin, ycos, x); break;
            case 6: _mp_real_sin_cos_opt_6(ysin, ycos, x); break;
            case 7: _mp_real_sin_cos_opt_7(ysin, ycos, x); break;
            case 8: _mp_real_sin_cos_opt_8(ysin, ycos, x); break;
            case 9: _mp_real_sin_cos_opt_9(ysin, ycos, x); break;
            default: _mp_real_sin_cos_opt_10(ysin, ycos, x); break;
        }
        *err = 6 * (ulong) _mp_real_sin_cos_default_r_inline(n) + 128;
        return 1;
    }

    if (n > MP_REAL_SIN_COS_OPT_MAX && n <= MP_REAL_DENSE_RT_MAX)
    {
        _mp_real_sin_cos_dense_rt(ysin, ycos, x, n);
        *err = 6 * (ulong) MP_REAL_DENSE_RT_R + 128;
        return 1;
    }

    return 0;
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

/* The exponential's reduction modulo log 2 (exp.c, and the nfloat
   exponentials): with X holding |x| < 2^(FLINT_BITS - 2) in a frame of
   N fraction limbs and one integral limb, q = _mp_real_exp_quotient of
   its top two limbs x1 + x0/B is |x| / log 2 or one less: with
   1/log 2 = 1 + R/B^2, R = (R1, R0) the top fraction limbs of 1/log 2,
   the product is x1 + h1 + (x0 + l1 + hi(x1 R0))/B plus positive terms
   below 2/B (h1 B + l1 = x1 R1), which are dropped before flooring. */
FLINT_FORCE_INLINE ulong
_mp_real_exp_quotient(ulong x1, ulong x0)
{
    nn_srcptr R = _mp_real_const_inv_log2_frac_static + MP_REAL_CONST_STATIC_N - 2;
    ulong h, l, c, u, s1, s0;

    umul_ppmm(h, l, x1, R[1]);
    umul_ppmm(c, u, x1, R[0]);
    (void) u;
    add_ssaaaa(s1, s0, UWORD(0), x0, UWORD(0), l);
    add_ssaaaa(s1, s0, s1, s0, UWORD(0), c);
    (void) s0;
    return x1 + h + s1;
}

/* The reduction step: X holds |x| (neg = 0) or its two's complement
   B^(N+1) - |x| (neg = 1) in the frame above, L = floor(log 2 B^N) and q
   from _mp_real_exp_quotient. One pass X - q L resp. X + (q + 1) L, and
   one more L if q was one low, leaves t = |x| - q L resp. (q + 1) L - |x|
   in [0, L] in X (X[N] = 0), and returns k with x = k log 2 + t: k = q
   resp. -(q + 1) (up to that correction). */
FLINT_FORCE_INLINE slong
_mp_real_exp_reduce_step(nn_ptr X, nn_srcptr L, slong N, ulong q, int neg)
{
    ulong cy;

    if (!neg)
    {
        MP_REAL_SUBMUL_1(cy, X, L, N, q);
        X[N] -= cy;
        if (X[N] != 0 || mpn_cmp(X, L, N) >= 0)
        {
            /* q was one low */
            X[N] -= mpn_sub_n(X, X, L, N);
            q++;
        }
        FLINT_ASSERT(X[N] == 0);
        return (slong) q;
    }
    else
    {
        q++;
        MP_REAL_ADDMUL_1(cy, X, L, N, q);
        X[N] += cy;
        if (X[N] != 0)
        {
            /* negative (q was one low): one more L */
            X[N] += mpn_add_n(X, X, L, N);
            q++;
        }
        FLINT_ASSERT(X[N] == 0);
        return -(slong) q;
    }
}

/* (y, n + 1) = exp(v) for v in [0, 1) at n fraction limbs: the kernel
   mp_real_exp_bits uses (the reduced series for small v, per-size,
   bitwise, diophantine or notab by n); err <= 9r + 100 with r <= 768
   the bitwise parameter, less for the others (exp.c) */
void _mp_real_exp_kernel(nn_ptr y, ulong * err, nn_srcptr v, slong n);

/* (y, n + 1) = exp(-v) for v in [0, 1) at n fraction limbs by the
   alternating series (or cosh - sinh) when v has enough leading zero
   bits for the series to beat the reduction and the kernel, returning
   1 (err <= 10); else returns 0 and leaves y alone (exp.c) */
int _mp_real_exp_neg_series(nn_ptr y, ulong * err, nn_srcptr v, slong n);

/* (y, n) = atan(v) for v in [0, 1) at n fraction limbs: the kernel
   mp_real_atan_bits uses (series for small v, per-size or bitwise, or
   Newton beyond 600 limbs); err <= 4r + 64 (atan.c) */
void _mp_real_atan_kernel(nn_ptr y, ulong * err, nn_srcptr v, slong n);

/* (y, n) = log1p(v) for v in [0, 1) at n <= 600 fraction limbs: the
   per-size or bitwise kernel of mp_real_log_bits; err <= 3r + 64
   (log.c) */
void _mp_real_log1p_kernel(nn_ptr y, ulong * err, nn_srcptr v, slong n);

/* the number of leading zero bits of s = |u - 1| from which
   log(1 +- s) = +-2 atanh(s / (2 +- s)) by the atanh series beats the
   log1p kernel, at n kernel limbs, for u above 1 (above = 1) resp.
   below (log.c) */
slong _mp_real_log_series_min_z(slong n, int above);

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

/* The exponent e2 of a remainder bound 2^e2 for a result of magnitude
   below 2^emid at prec bits, raised to emid - prec - 2 FLINT_BITS when
   smaller: still a valid bound, of negligible relative size, which pads
   the mantissa by at most a few limbs beyond prec bits rather than down
   to 2^e2 (B^(2^55) limbs for the tiniest arguments). */
FLINT_FORCE_INLINE slong
_mp_real_err_exp_clamp(slong e2, slong emid, slong prec)
{
    return FLINT_MAX(e2, emid - prec - 2 * FLINT_BITS);
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

/* res = w sum_{k=1}^{N} s_k c_{k-1} v^(k-1) to p limbs, v = w (odd = 0) or
   w^2 (odd = 1), s_k = (-1)^(k+1) if alternating else 1, |w| < 2^uexp,
   N <= 16 (newton.c); the caller divides by the common denominator of
   the coefficients and adds the tail */
void _mp_real_newton_series(mp_real_t res, const mp_real_t w, slong uexp,
    slong p, const ulong * c, slong N, int odd, int alternating);

#ifdef __cplusplus
}
#endif

#endif
