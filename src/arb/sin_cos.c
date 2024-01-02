/*
    Copyright (C) 2012-2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "arb.h"

#define TMP_ALLOC_LIMBS(__n) TMP_ALLOC((__n) * sizeof(mp_limb_t))
#define MAGLIM(prec) FLINT_MAX(65536, (4*prec))

static void
mag_nonzero_fast_mul(mag_t z, const mag_t x, const mag_t y)
{
    MAG_MAN(z) = MAG_FIXMUL(MAG_MAN(x), MAG_MAN(y)) + LIMB_ONE;
    MAG_EXP(z) = MAG_EXP(x) + MAG_EXP(y);
    MAG_FAST_ADJUST_ONE_TOO_SMALL(z);
}

static void
mag_nonzero_fast_add(mag_t z, const mag_t x, const mag_t y)
{
    slong shift = MAG_EXP(x) - MAG_EXP(y);

    if (shift == 0)
    {
        MAG_EXP(z) = MAG_EXP(x);
        MAG_MAN(z) = MAG_MAN(x) + MAG_MAN(y);
        MAG_FAST_ADJUST_ONE_TOO_LARGE(z); /* may need two adjustments */
    }
    else if (shift > 0)
    {
        MAG_EXP(z) = MAG_EXP(x);

        if (shift >= MAG_BITS)
            MAG_MAN(z) = MAG_MAN(x) + LIMB_ONE;
        else
            MAG_MAN(z) = MAG_MAN(x) + (MAG_MAN(y) >> shift) + LIMB_ONE;
    }
    else
    {
        shift = -shift;
        MAG_EXP(z) = MAG_EXP(y);

        if (shift >= MAG_BITS)
            MAG_MAN(z) = MAG_MAN(y) + LIMB_ONE;
        else
            MAG_MAN(z) = MAG_MAN(y) + (MAG_MAN(x) >> shift) + LIMB_ONE;
    }

    MAG_FAST_ADJUST_ONE_TOO_LARGE(z);
}

static int
mag_nonzero_fast_cmp(const mag_t x, const mag_t y)
{
    if (MAG_EXP(x) == MAG_EXP(y))
        return (MAG_MAN(x) < MAG_MAN(y)) ? -1 : 1;
    else
        return (MAG_EXP(x) < MAG_EXP(y)) ? -1 : 1;
}

static void
mag_fast_set(mag_t x, const mag_t y)
{
    MAG_EXP(x) = MAG_EXP(y);
    MAG_MAN(x) = MAG_MAN(y);
}

void
_arb_sin_cos(arb_t zsin, arb_t zcos, const arf_t x, const mag_t xrad, slong prec)
{
    int want_sin, want_cos;
    slong radexp, exp, wp, wn, N, r, wprounded, maglim, orig_prec;
    mp_ptr tmp, w, sina, cosa, sinb, cosb, ta, tb;
    mp_ptr sinptr, cosptr;
    mp_limb_t p1, q1bits, p2, q2bits, error, error2, p1_tab1, radman;
    int negative, inexact, octant;
    int sinnegative, cosnegative, swapsincos;
    TMP_INIT;

    /* PART 1: special cases and setup. */
    orig_prec = prec;

    /* Below, both x and xrad will be finite, and x will be nonzero. */
    if (mag_is_inf(xrad) || arf_is_special(x))
    {
        _arb_sin_cos_generic(zsin, zcos, x, xrad, prec);
        return;
    }

    exp = ARF_EXP(x);
    maglim = MAGLIM(prec);
    negative = ARF_SGNBIT(x);

    /* Unlikely: tiny or huge midpoint. As written, this test also
       catches any bignum exponents. */
    if (exp < -(prec/2) - 2 || exp > maglim)
    {
        _arb_sin_cos_generic(zsin, zcos, x, xrad, prec);
        return;
    }

    want_sin = (zsin != NULL);
    want_cos = (zcos != NULL);

    /* Copy the radius data. */
    radexp = MAG_EXP(xrad);
    radman = MAG_MAN(xrad);

    if (radman != 0)
    {
        /* Clamp the radius exponent to a safe range. */
        if (radexp < MAG_MIN_LAGOM_EXP || radexp > MAG_MAX_LAGOM_EXP)
        {
            /* Very wide... */
            if (fmpz_sgn(MAG_EXPREF(xrad)) > 0)
            {
                _arb_sin_cos_wide(zsin, zcos, x, xrad, prec);
                return;
            }

            radman = MAG_ONE_HALF;
            radexp = MAG_MIN_LAGOM_EXP + 1;
        }

        /* Use wide algorithm. */
        if (radexp >= -24)
        {
            _arb_sin_cos_wide(zsin, zcos, x, xrad, prec);
            return;
        }

        /* Regular case: decrease precision to match generic max. accuracy. */
        /* Note: near x=0, the error can be quadratic for cos. */
        if (want_cos && exp < -2)
            prec = FLINT_MIN(prec, 20 - FLINT_MAX(exp, radexp) - radexp);
        else
            prec = FLINT_MIN(prec, 20 - radexp);
    }

    if (prec >= ARB_SIN_COS_ATAN_REDUCTION_PREC && prec <= ARB_SIN_COS_ATAN_REDUCTION_DEFAULT_MAX_PREC)
    {
        mag_t err;
        mag_init(err);
        mag_set(err, xrad);
        arb_sin_cos_arf_atan_reduction(zsin, zcos, x, prec);
        if (zsin != NULL) arb_add_error_mag(zsin, err);
        if (zcos != NULL) arb_add_error_mag(zcos, err);
        mag_clear(err);
        return;
    }

    /* Absolute working precision (NOT rounded to a limb multiple) */
    wp = prec + 8;
    if (want_sin && exp <= 0)
        wp += (-exp);
    /* Number of limbs */
    wn = (wp + FLINT_BITS - 1) / FLINT_BITS;
    /* Precision rounded to a number of bits */
    wprounded = FLINT_BITS * wn;
    /* Don't be close to the boundary (to allow adding adding the
       Taylor series truncation error without overflow) */
    wp = FLINT_MAX(wp, wprounded - (FLINT_BITS - 4));

    /* Too high precision to use table -- use generic algorithm */
    if (wp > ARB_SIN_COS_TAB2_PREC)
    {
        _arb_sin_cos_generic(zsin, zcos, x, xrad, orig_prec);
        return;
    }

    /* PART 2: the actual computation. */

    TMP_START;
    tmp = TMP_ALLOC_LIMBS(9 * wn);
    w    = tmp;         /* requires wn limbs */
    sina = w    + wn;   /* requires wn limbs */
    cosa = sina + wn;   /* requires wn limbs */
    sinb = cosa + wn;   /* requires wn limbs */
    cosb = sinb + wn;   /* requires wn limbs */
    ta   = cosb + wn;   /* requires 2*wn limbs */
    tb   = ta + 2*wn;   /* requires 2*wn limbs */

    /* reduce modulo pi/4 */
    if (_arb_get_mpn_fixed_mod_pi4(w, NULL, &octant, &error, x, wn) == 0)
    {
        /* may run out of precision for pi/4 */
        _arb_sin_cos_generic(zsin, zcos, x, xrad, orig_prec);
        TMP_END;
        return;
    }

    sinnegative = (octant >= 4) ^ negative;
    cosnegative = (octant >= 2 && octant <= 5);
    swapsincos = (octant == 1 || octant == 2 || octant == 5 || octant == 6);

    /* Table-based argument reduction (1 or 2 steps) */
    if (wp <= ARB_SIN_COS_TAB1_PREC)
    {
        q1bits = ARB_SIN_COS_TAB1_BITS;
        q2bits = 0;

        p1 = p1_tab1 = w[wn-1] >> (FLINT_BITS - q1bits);
        w[wn-1] -= (p1 << (FLINT_BITS - q1bits));
        p2 = 0;

        /* p1_tab1 will be used for the error bounds at the end. */
        p1_tab1 = p1;
    }
    else
    {
        q1bits = ARB_SIN_COS_TAB21_BITS;
        q2bits = ARB_SIN_COS_TAB21_BITS + ARB_SIN_COS_TAB22_BITS;

        /* p1_tab1 will be used for the error bounds at the end. */
        p1_tab1 = w[wn-1] >> (FLINT_BITS - ARB_SIN_COS_TAB1_BITS);

        p1 = w[wn-1] >> (FLINT_BITS - q1bits);
        w[wn-1] -= (p1 << (FLINT_BITS - q1bits));
        p2 = w[wn-1] >> (FLINT_BITS - q2bits);
        w[wn-1] -= (p2 << (FLINT_BITS - q2bits));
    }

    /* |w| <= 2^-r */
    r = _arb_mpn_leading_zeros(w, wn);

    /* Choose number of terms N such that Taylor series truncation
       error is <= 2^-wp */
    N = _arb_exp_taylor_bound(-r, wp);

    /* the summation for sin/cos is actually done to (2N-1)! */
    N = (N + 1) / 2;

    if (N < 14)
    {
        /* Evaluate Taylor series */
        _arb_sin_cos_taylor_rs(sina, cosa, &error2, w, wn, N, 0, 1);
        /* Taylor series evaluation error */
        error += error2;
        /* Taylor series truncation error */
        error += UWORD(1) << (wprounded-wp);
    }
    else  /* Compute cos(a) from sin(a) using a square root. */
    {
        /* Evaluate Taylor series */
        _arb_sin_cos_taylor_rs(sina, cosa, &error2, w, wn, N, 1, 1);
        error += error2;
        error += UWORD(1) << (wprounded-wp);

        if (flint_mpn_zero_p(sina, wn))
        {
            flint_mpn_store(cosa, wn, LIMB_ONES);
            error = FLINT_MAX(error, 1);
        }
        else
        {
            flint_mpn_sqr(ta, sina, wn);
            /* 1 - s^2 (negation guaranteed to have borrow) */
            mpn_neg(ta, ta, 2 * wn);
            /* top limb of ta must be nonzero, so no need to normalize
               before calling mpn_sqrtrem */
            mpn_sqrtrem(cosa, ta, ta, 2 * wn);

            /* When converting sin to cos, the error for cos must be
               smaller than the error for sin; but we also get 1 ulp
               extra error from the square root. */
            error += 1;
        }
    }

    /*
    sin(a+b) = sin(a)*cos(b) + cos(a)*sin(b)
    cos(a+b) = cos(a)*cos(b) - sin(a)*sin(b)

    [F1+e1]*[F2+e2] + [F3+e3]*[F4+e4] - F1 F2 + F3 F4 =

       e1 e2 + e3 e4 + e2 F1 + e1 F2 + e4 F3 + e3 F4

       <= (e1 + e2 + e3 + e4) + 1    (ulp)

       <= 2*left_err + 2*right_err + 1     (ulp)

    Truncating both terms before adding adds another 2 ulp, so the error
    is bounded by

       <= 2*left_err + 2*right_err + 3     (ulp)
    */

    if (p1 == 0 && p2 == 0)  /* no table lookups */
    {
        sinptr = sina;
        cosptr = cosa;
    }
    else if (p1 == 0 || p2 == 0)    /* only one table lookup */
    {
        mp_srcptr sinc, cosc;

        if (wp <= ARB_SIN_COS_TAB1_PREC)  /* must be in table 1 */
        {
            sinc = arb_sin_cos_tab1[2 * p1] + ARB_SIN_COS_TAB1_LIMBS - wn;
            cosc = arb_sin_cos_tab1[2 * p1 + 1] + ARB_SIN_COS_TAB1_LIMBS - wn;
        }
        else if (p1 != 0)
        {
            sinc = arb_sin_cos_tab21[2 * p1] + ARB_SIN_COS_TAB2_LIMBS - wn;
            cosc = arb_sin_cos_tab21[2 * p1 + 1] + ARB_SIN_COS_TAB2_LIMBS - wn;
        }
        else
        {
            sinc = arb_sin_cos_tab22[2 * p2] + ARB_SIN_COS_TAB2_LIMBS - wn;
            cosc = arb_sin_cos_tab22[2 * p2 + 1] + ARB_SIN_COS_TAB2_LIMBS - wn;
        }

        if ((want_sin && !swapsincos) || (want_cos && swapsincos))
        {
            flint_mpn_mul_n(ta, sina, cosc, wn);
            flint_mpn_mul_n(tb, cosa, sinc, wn);
            mpn_add_n(w, ta + wn, tb + wn, wn);
        }

        if ((want_cos && !swapsincos) || (want_sin && swapsincos))
        {
            flint_mpn_mul_n(ta, cosa, cosc, wn);
            flint_mpn_mul_n(tb, sina, sinc, wn);
            mpn_sub_n(ta, ta + wn, tb + wn, wn);
        }

        sinptr = w;
        cosptr = ta;

        error = 2 * error + 2 * 1 + 3;
    }
    else        /* two table lookups, must be in table 2 */
    {
        mp_srcptr sinc, cosc, sind, cosd;

        sinc = arb_sin_cos_tab21[2 * p1] + ARB_SIN_COS_TAB2_LIMBS - wn;
        cosc = arb_sin_cos_tab21[2 * p1 + 1] + ARB_SIN_COS_TAB2_LIMBS - wn;
        sind = arb_sin_cos_tab22[2 * p2] + ARB_SIN_COS_TAB2_LIMBS - wn;
        cosd = arb_sin_cos_tab22[2 * p2 + 1] + ARB_SIN_COS_TAB2_LIMBS - wn;

        flint_mpn_mul_n(ta, sinc, cosd, wn);
        flint_mpn_mul_n(tb, cosc, sind, wn);
        mpn_add_n(sinb, ta + wn, tb + wn, wn);

        flint_mpn_mul_n(ta, cosc, cosd, wn);
        flint_mpn_mul_n(tb, sinc, sind, wn);
        mpn_sub_n(cosb, ta + wn, tb + wn, wn);

        error2 = 2 * 1 + 2 * 1 + 3;

        if ((want_sin && !swapsincos) || (want_cos && swapsincos))
        {
            flint_mpn_mul_n(ta, sina, cosb, wn);
            flint_mpn_mul_n(tb, cosa, sinb, wn);
            mpn_add_n(w, ta + wn, tb + wn, wn);
        }

        if ((want_cos && !swapsincos) || (want_sin && swapsincos))
        {
            flint_mpn_mul_n(ta, cosa, cosb, wn);
            flint_mpn_mul_n(tb, sina, sinb, wn);
            mpn_sub_n(ta, ta + wn, tb + wn, wn);
        }

        error = 2 * error + 2 * error2 + 3;

        sinptr = w;
        cosptr = ta;
    }

    /* PART 3: compute propagated error and write output */

    if (swapsincos)
    {
        mp_ptr tmptr = sinptr;
        sinptr = cosptr;
        cosptr = tmptr;
    }

    /*
    We have two sources of error.

    1. Computation error:
        error * 2^(-wprounded)

    2. With input radius r != 0, the propagated error bound:
        sin(x):  min(2, r, |cos(x)|*r  +  0.5*r^2)
        cos(x):  min(2, r, |sin(x)|*r  +  0.5*r^2)

       We skip the min by 2 since this is unnecessary with the
       separate code for wide intervals.
    */
    if (radman == 0)
    {
        if (want_sin)
        {
            mag_set_ui_2exp_si(arb_radref(zsin), error, -wprounded);
            if (want_cos)
                mag_set(arb_radref(zcos), arb_radref(zsin));
        }
        else
        {
            mag_set_ui_2exp_si(arb_radref(zcos), error, -wprounded);
        }
    }
    else
    {
        mag_t sin_err, cos_err, quadratic, comp_err, xrad_copy;
        mp_limb_t A_sin, A_cos, A_exp;

        /* Copy xrad to support aliasing (note: the exponent has
           also been clamped earlier). */
        MAG_MAN(xrad_copy) = radman;
        MAG_EXP(xrad_copy) = radexp;

        /* Bound computed error. */
        if (error != 0)
        {
            mag_init(comp_err); /* no need to free */
            mag_set_ui_2exp_si(comp_err, error, -wprounded);
        }

        /* Bound quadratic term for propagated error: 0.5*r^2 */
        mag_init(quadratic); /* no need to free */
        mag_nonzero_fast_mul(quadratic, xrad_copy, xrad_copy);
        MAG_EXP(quadratic) -= 1;

        /* Bound linear term for propagated error: cos(x)*r, sin(x)*r. */
        /* Note: we could have used the computed values, but then we would
           need to incorporate the computed error which would be slightly
           messier, and we would also need extra cases when only computing
           one of the functions. */
        A_cos = arb_sin_cos_tab1[2 * p1_tab1][ARB_SIN_COS_TAB1_LIMBS - 1];
        A_sin = arb_sin_cos_tab1[2 * p1_tab1 + 1][ARB_SIN_COS_TAB1_LIMBS - 1];

        /* Note: ARB_SIN_COS_TAB1_BITS == 8 */
        /* Adding 2 ulps (here ulp = 2^-8) gives an upper bound.
           The truncated table entry underestimates the sine or
           cosine of x by at most 1 ulp, and the top bits of x
           underestimate x by at most 1 ulp. */
        A_sin = (A_sin >> (FLINT_BITS - ARB_SIN_COS_TAB1_BITS)) + 2;
        A_cos = (A_cos >> (FLINT_BITS - ARB_SIN_COS_TAB1_BITS)) + 2;
        A_exp = -ARB_SIN_COS_TAB1_BITS;
        if (swapsincos)
        {
            mp_limb_t tt = A_sin;
            A_sin = A_cos;
            A_cos = tt;
        }
        A_sin *= ((MAG_MAN(xrad_copy) >> (MAG_BITS - ARB_SIN_COS_TAB1_BITS)) + LIMB_ONE);
        A_cos *= ((MAG_MAN(xrad_copy) >> (MAG_BITS - ARB_SIN_COS_TAB1_BITS)) + LIMB_ONE);
        A_exp -= ARB_SIN_COS_TAB1_BITS;
        A_exp += radexp;

        if (want_sin)
        {
            mag_init(sin_err);
            mag_set_ui_2exp_si(sin_err, A_sin, A_exp);
            mag_nonzero_fast_add(sin_err, sin_err, quadratic);

            /* The propagated error is certainly at most r */
            if (mag_nonzero_fast_cmp(sin_err, xrad_copy) > 0)
                mag_fast_set(sin_err, xrad_copy);

            /* Add the computed error. */
            if (error != 0)
                mag_nonzero_fast_add(sin_err, sin_err, comp_err);

            /* Set it, and clear the original output variable which could
               have a bignum exponent. */
            mag_swap(arb_radref(zsin), sin_err);
            mag_clear(sin_err);
        }

        /* The same as above. */
        if (want_cos)
        {
            mag_init(cos_err);
            mag_set_ui_2exp_si(cos_err, A_cos, A_exp);
            mag_nonzero_fast_add(cos_err, cos_err, quadratic);

            if (mag_nonzero_fast_cmp(cos_err, xrad_copy) > 0)
                mag_fast_set(cos_err, xrad_copy);

            if (error != 0)
                mag_nonzero_fast_add(cos_err, cos_err, comp_err);
            mag_swap(arb_radref(zcos), cos_err);
            mag_clear(cos_err);
        }
    }

    /* Set the midpoints. */
    if (want_sin)
    {
        inexact = _arf_set_mpn_fixed(arb_midref(zsin), sinptr,
            wn, wn, sinnegative, prec, ARB_RND);
        if (inexact)
            arf_mag_add_ulp(arb_radref(zsin),
                arb_radref(zsin), arb_midref(zsin), prec);
    }

    if (want_cos)
    {
        inexact = _arf_set_mpn_fixed(arb_midref(zcos), cosptr,
            wn, wn, cosnegative, prec, ARB_RND);
        if (inexact)
            arf_mag_add_ulp(arb_radref(zcos),
                arb_radref(zcos), arb_midref(zcos), prec);
    }

    TMP_END;
}

void
arb_sin_cos(arb_t s, arb_t c, const arb_t x, slong prec)
{
    _arb_sin_cos(s, c, arb_midref(x), arb_radref(x), prec);
}

void
arb_sin(arb_t s, const arb_t x, slong prec)
{
    _arb_sin_cos(s, NULL, arb_midref(x), arb_radref(x), prec);
}

void
arb_cos(arb_t c, const arb_t x, slong prec)
{
    _arb_sin_cos(NULL, c, arb_midref(x), arb_radref(x), prec);
}
