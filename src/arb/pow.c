/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb/elem.h"

#define BINEXP_LIMIT 64

typedef enum {POSITIVE = 0, NEGATIVE_EVEN, NEGATIVE_ODD} sign_type;

/* |x|^y for exact x, y by the kernel _mp_real_pow_fp (exp(y log |x|) in
   one evaluation): 1 if done, 0 when the arguments are outside the
   kernel's range (the exponents must be words; |y log |x|| below
   2^maglim, where exp gives up: [0 +- inf] resp. a tiny ball) */
static int
_arb_pow_arf(arb_t z, const arf_t x, const arf_t y, slong prec)
{
    nn_srcptr xd, yd;
    slong xn, yn, xe, ye, nk, qn, maglim = FLINT_MAX(128, 2 * prec);
    nn_ptr out, q;
    ulong err;
    int code, tsgn;
    TMP_INIT;

    if (arf_is_special(x) || arf_is_special(y) || COEFF_IS_MPZ(ARF_EXP(x)) || COEFF_IS_MPZ(ARF_EXP(y)))
        return 0;

    xe = ARF_EXP(x);
    ye = ARF_EXP(y);
    if (xe >= (WORD(1) << (FLINT_BITS - 3)) || xe <= -(WORD(1) << (FLINT_BITS - 3)))
        return 0;

    ARF_GET_MPN_READONLY(xd, xn, x);
    ARF_GET_MPN_READONLY(yd, yn, y);

    /* |x| = 1 */
    if (xn == 1 && xe == 1 && xd[0] == LIMB_TOP)
    {
        arb_one(z);
        return 1;
    }

#if FLINT_BITS == 64
    /* the register paths at one and two limbs (arf/elem.h), shorter
       mantissas padded */
    if (prec <= 2 * FLINT_BITS && xn <= 2 && yn <= 2)
    {
        ulong y3 = 0, y2 = 0, y1 = 0, y0 = 0;
        slong k = 0;

        if (xn == 1 && yn == 1 && prec <= FLINT_BITS)
        {
            if (_arf_small_pow_1(&y2, &y1, &y0, &k, &err, xd[0], xe, yd[0], ye, ARF_SGNBIT(y)) != ARF_ELEM_UNCERTAIN)
            {
                _arb_small_set_1(z, y2, y1, y0, k, err, 0, prec);
                return 1;
            }
        }
        else
        {
            if (_arf_small_pow_2(&y3, &y2, &y1, &y0, &k, &err, xd[xn - 1], (xn == 2) ? xd[0] : 0, xe,
                    yd[yn - 1], (yn == 2) ? yd[0] : 0, ye, ARF_SGNBIT(y)) != ARF_ELEM_UNCERTAIN)
            {
                _arb_small_set_2(z, y3, y2, y1, y0, k, err, 0, prec);
                return 1;
            }
        }
    }
#endif

    nk = _arb_elem_limbs(prec);
    TMP_START;
    out = TMP_ALLOC((nk + 2 + (maglim + FLINT_BITS - 1) / FLINT_BITS + 2) * sizeof(ulong));
    q = out + nk + 2;

    code = _mp_real_pow_fp(out, &err, q, &qn, &tsgn, xd, xn, xe, yd, yn, ye, ARF_SGNBIT(y), _arb_elem_bits(prec), maglim);

    if (code == MP_REAL_POW_HUGE)
    {
        if (!tsgn)
            arb_zero_pm_inf(z);
        else
        {
            /* exp(t) < 2^(-2^maglim) */
            fmpz_t t;
            fmpz_init(t);
            fmpz_set_si(t, -1);
            fmpz_mul_2exp(t, t, maglim);
            arf_one(arb_midref(z));
            mag_one(arb_radref(z));
            arb_mul_2exp_fmpz(z, z, t);
            fmpz_clear(t);
        }
    }
    else if (code == MP_REAL_POW_TINY)
    {
        /* 1 + delta, |delta| < 2^-(FLINT_BITS nk + 3) */
        arf_one(arb_midref(z));
        mag_set_ui_2exp_si(arb_radref(z), 1, -FLINT_BITS * nk - 3);
    }
    else if (code == MP_REAL_POW_ONE)
        arb_one(z);
    else
    {
        _arb_set_mpn_err(z, out, nk + 1, -FLINT_BITS * nk, err, err, 0, prec);
        if (qn != 0)
        {
            fmpz_t qz;
            fmpz_init(qz);
            fmpz_set_ui_array(qz, q, qn);
            if (code)
                fmpz_neg(qz, qz);
            arb_mul_2exp_fmpz(z, z, qz);
            fmpz_clear(qz);
        }
    }

    TMP_END;
    return 1;
}

static void
_arb_pow_exp(arb_t z, const arb_t x, sign_type negx, const arb_t y, slong prec)
{
    arb_t t;

    if (arb_is_exact(x) && arb_is_exact(y) && (negx != POSITIVE || arf_sgn(arb_midref(x)) > 0))
    {
        arf_t ax;
        int ok;
        arf_init(ax);
        arf_abs(ax, arb_midref(x));
        ok = _arb_pow_arf(z, ax, arb_midref(y), prec);
        arf_clear(ax);
        if (ok)
        {
            if (negx == NEGATIVE_ODD)
                arb_neg(z, z);
            return;
        }
    }

    arb_init(t);

    if (negx == POSITIVE)
        arb_log(t, x, prec);
    else
    {
        arb_neg(t, x);
        arb_log(t, t, prec);
    }

    arb_mul(t, t, y, prec);
    arb_exp(z, t, prec);

    if (negx == NEGATIVE_ODD)
        arb_neg(z, z);

    arb_clear(t);
}

/* x^n over a wide ball x (a relative radius above 1/(4n), or 0 inside)
   for the exact integer n != 0, as the hull of the range: x^n is
   monotone on each side of 0 (increasing for odd n, even powers
   increasing in |x|), so the range is spanned by the endpoints, or
   by 0 and the larger endpoint in absolute value.  Powers of the
   endpoints by arb_pow_fmpz_binexp (exact arguments: tight) up to
   2^BINEXP_LIMIT, beyond as exp(n log) of the exact endpoints.  A ball
   containing 0 with n < 0 gives an indeterminate result. */
static void
_arb_pow_fmpz_wide(arb_t z, const arb_t x, const fmpz_t n, slong prec)
{
    arf_t a, b;
    int odd = fmpz_is_odd(n), zero_inside = arb_contains_zero(x);
    int neg = (fmpz_sgn(n) < 0);

    if (zero_inside && neg)
    {
        arb_indeterminate(z);
        return;
    }

    arf_init(a);
    arf_init(b);
    arb_get_interval_arf(a, b, x, prec + 10);

    if (!odd)
    {
        /* even: the range of |x|^n; |x| over [0 or min(|a|, |b|), max(|a|, |b|)] */
        arf_abs(a, a);
        arf_abs(b, b);
        if (arf_cmp(a, b) > 0)
            arf_swap(a, b);
        if (zero_inside)
            arf_zero(a);
    }

    {
        /* the powers of the (exact) endpoints as balls, and their hull:
           binary exponentiation for small n, else exp(n log) through
           arb_pow itself (on exact arguments) */
        arb_t s, t, e;
        arb_init(s);
        arb_init(t);
        arb_set_arf(s, a);
        arb_set_arf(t, b);
        if (fmpz_bits(n) <= BINEXP_LIMIT)
        {
            arb_pow_fmpz_binexp(s, s, n, prec + 10);
            arb_pow_fmpz_binexp(t, t, n, prec + 10);
        }
        else
        {
            arb_init(e);
            arb_set_fmpz(e, n);
            arb_pow(s, s, e, prec + 10);
            arb_pow(t, t, e, prec + 10);
            arb_clear(e);
        }
        if (neg)
            arb_swap(s, t);
        arb_union(z, s, t, prec);
        arb_clear(s);
        arb_clear(t);
    }

    arf_clear(a);
    arf_clear(b);
}

/* whether the ball x = [m +- r] is wide for the integer power n: r >= |m|
   / (4 |n|) (the hull of the range is then clearly tighter than the ball
   arithmetic of binary exponentiation), or 0 inside */
static int
_arb_pow_is_wide(const arb_t x, const fmpz_t n)
{
    mag_t t;
    int res;

    if (mag_is_zero(arb_radref(x)))
        return 0;
    if (arb_contains_zero(x))
        return 1;

    mag_init(t);
    arf_get_mag_lower(t, arb_midref(x));
    if (fmpz_bits(n) <= 60)
        mag_div_ui(t, t, 4 * (ulong) FLINT_ABS(fmpz_get_si(n)));
    else
        mag_zero(t);
    res = (mag_cmp(arb_radref(x), t) >= 0);
    mag_clear(t);
    return res;
}

/* whether the ball y may contain an integer */
static int
_arb_contains_int(const arb_t y)
{
    arf_t lo, hi;
    int res;
    if (!arb_is_finite(y))
        return 1;
    arf_init(lo);
    arf_init(hi);
    arb_get_lbound_arf(lo, y, MAG_BITS);
    arb_get_ubound_arf(hi, y, MAG_BITS);
    /* an integer in [lo, hi] iff floor(hi) >= ceil(lo) */
    arf_floor(hi, hi);
    arf_ceil(lo, lo);
    res = arf_cmp(hi, lo) >= 0;
    arf_clear(lo);
    arf_clear(hi);
    return res;
}

void
arb_pow(arb_t z, const arb_t x, const arb_t y, slong prec)
{
    sign_type s;

    if (arb_is_zero(y))
    {
        arb_one(z);
        return;
    }

    if (arb_is_zero(x))
    {
        if (arb_is_positive(y))
            arb_zero(z);
        else
            arb_indeterminate(z);
        return;
    }

    s = POSITIVE;

    if (arb_is_exact(y) && arf_is_finite(arb_midref(x)))
    {
        const arf_struct * ymid = arb_midref(y);

        if (arf_is_int(ymid) && arf_cmpabs_2exp_si(ymid, 1024) < 0)
        {
            /* an integer: binary exponentiation of a precise ball (small
               exponents), the hull of the range over a wide one */
            fmpz_t e;
            fmpz_init(e);
            arf_get_fmpz(e, ymid, ARF_RND_DOWN);

            if (_arb_pow_is_wide(x, e))
                _arb_pow_fmpz_wide(z, x, e, prec);
            else if (arf_cmpabs_2exp_si(ymid, BINEXP_LIMIT) < 0)
                arb_pow_fmpz_binexp(z, x, e, prec);
            else if (arf_sgn(arb_midref(x)) < 0)
            {
                /* (-x)^n = (-1)^n x^n (exp(n log(-x)) below) */
                s = fmpz_is_even(e) ? NEGATIVE_EVEN : NEGATIVE_ODD;
                fmpz_clear(e);
                goto generic;
            }
            else
            {
                fmpz_clear(e);
                goto generic;
            }

            fmpz_clear(e);
            return;
        }
        else if (arf_cmpabs_2exp_si(ymid, BINEXP_LIMIT) < 0 &&
            arf_is_int_2exp_si(ymid, -1) && !arb_contains_negative(x))
        {
            /* a small half-integer: (sqrt x)^(2y) (a base reaching into
               the negative numbers goes below, to the range over its
               nonnegative part) */
            fmpz_t e;
            fmpz_init(e);
            arf_get_fmpz_fixed_si(e, ymid, -1);
            if (fmpz_sgn(e) >= 0)
            {
                arb_sqrt(z, x, prec + fmpz_bits(e));
                arb_pow_fmpz_binexp(z, z, e, prec);
            }
            else
            {
                fmpz_neg(e, e);
                arb_rsqrt(z, x, prec + fmpz_bits(e));
                arb_pow_fmpz_binexp(z, z, e, prec);
            }
            fmpz_clear(e);
            return;
        }
        else if (arf_is_int(ymid) && arf_sgn(arb_midref(x)) < 0)
            /* use (-x)^n = (-1)^n * x^n to avoid NaNs
               at least at high enough precision */
            s = arf_is_int_2exp_si(ymid, 1) ? NEGATIVE_EVEN : NEGATIVE_ODD;
        /* Fallthrough */
    }

generic:
    /* (the negative part of x may be ignored only when y contains no
       integer: (-3)^1 = -3 is a real value of [-1 +- 2]^[1 +- 1/2]) */
    if (arb_contains_zero(x) && arb_is_nonnegative(y) && arb_contains_positive(x)
        && (!arb_contains_negative(x) || !_arb_contains_int(y)))
    {
        /* x contains 0 (and positive numbers), y is nonnegative: over the
           domain x >= 0, x^y is nondecreasing in x and 0 at x = 0 (for
           y > 0; x^0 = 1 lies below the bound at y = 0 too), so the range
           is [0, R^y] for the upper endpoint R of x (the negative part of
           x has no real powers).  The upper bound is the upper endpoint
           of t^y for t = [R/2 +- R/2] (which excludes 0 up to its endpoint
           and so takes the generic route), and the result [u/2 +- u/2];
           this costs the precision of the radius: MAG_BITS suffice. */
        arb_t t;
        arf_t R;

        prec = FLINT_MIN(prec, MAG_BITS);

        arf_init(R);
        arb_init(t);
        arb_get_ubound_arf(R, x, prec);
        arf_mul_2exp_si(R, R, -1);
        arf_mul_ui(arb_midref(t), R, 3, prec, ARF_RND_UP);
        arf_mul_2exp_si(arb_midref(t), arb_midref(t), -1);
        arf_get_mag(arb_radref(t), R);
        mag_mul_2exp_si(arb_radref(t), arb_radref(t), -1);
        arf_clear(R);

        arb_pow(z, t, y, prec);
        arb_clear(t);

        /* the upper bound u of z (rounded up), the result [u/2 +- u/2] */
        arb_get_ubound_arf(arb_midref(z), z, prec);
        arf_mul_2exp_si(arb_midref(z), arb_midref(z), -1);
        arf_get_mag(arb_radref(z), arb_midref(z));
        arf_set_mag(arb_midref(z), arb_radref(z));
        return;
    }

    _arb_pow_exp(z, x, s, y, prec);
}
