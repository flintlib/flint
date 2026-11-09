/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb/elem.h"

/* tanh (arb/elem.h): the point evaluation _arb_sinh_cosh_arf; a narrow
   radius (below 2^-10) propagated with tanh' = 1 - tanh^2 = sech^2,
   decreasing in |x|: over [m +- r], at most sech^2(|m| - r) <=
   sech^2(m) e^(2r) <= 2 (1 - |tanh m|) (1 + 4r); a wide ball with
   |m| >= 1 and r < 1/2 through exp(-2|x|) in ball arithmetic, otherwise
   by the values at the endpoints (tanh is increasing), which is sharp */
/* |m| - r >= 1 for the ball x = [m +- r] */
static int
_arb_tanh_ball_beyond_one(const arb_t x)
{
    arf_t t;
    int res;
    arf_init(t);
    arf_set_mag(t, arb_radref(x));
    if (arf_sgn(arb_midref(x)) < 0)
    {
        arf_add(t, t, arb_midref(x), MAG_BITS, ARF_RND_CEIL);   /* r + m = r - |m| */
        arf_neg(t, t);
    }
    else
        arf_sub(t, arb_midref(x), t, MAG_BITS, ARF_RND_FLOOR);  /* |m| - r */
    res = arf_cmp_ui(t, 1) >= 0;    /* |m| - r >= 1 */
    arf_clear(t);
    return res;
}

void
arb_tanh(arb_t y, const arb_t x, slong prec)
{
    slong radexp, mexp;
    mag_t r, t, u;

    if (mag_is_zero(arb_radref(x)))
    {
        _arb_sinh_cosh_arf(NULL, NULL, y, NULL, arb_midref(x), prec, FLINT_MAX(128, 2 * prec));
        return;
    }

    if (!arb_is_finite(x))
    {
        if (arf_is_nan(arb_midref(x)))
            arb_indeterminate(y);
        else if (arf_is_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)))
            arb_set_si(y, arf_sgn(arb_midref(x)));
        else
            arb_zero_pm_one(y);
        return;
    }

    radexp = MAG_EXP(arb_radref(x));
    if (COEFF_IS_MPZ(radexp))
        radexp = (fmpz_sgn(MAG_EXPREF(arb_radref(x))) < 0) ? -prec - 100 : 10;

    mexp = ARF_EXP(arb_midref(x));
    if (COEFF_IS_MPZ(mexp))
        mexp = (fmpz_sgn(ARF_EXPREF(arb_midref(x))) < 0) ? -prec - 100 : 100;

    /* a ball around 0 reaching beyond +-prec: [0 +- 1] (the ends within
       2 e^(-prec) of +-1) */
    if (radexp >= mexp + 2 && radexp >= (slong) FLINT_BIT_COUNT(prec) + 2)
    {
        arb_zero_pm_one(y);
        return;
    }

    if (radexp >= -9 && mexp >= 1
        && ((mexp >= 3 && radexp <= mexp - 2) || _arb_tanh_ball_beyond_one(x)))
    {
        /* a wide ball within |x| >= 1: tanh |x| = 1 - 2 E / (1 + E)
           with E = exp(-2|x|) as a ball, whose relative accuracy (the
           width of the argument, which arb_exp lets bound its working
           precision) is all the result needs: nearly sharp, and cheap */
        arb_t t;
        int neg = arf_sgn(arb_midref(x)) < 0;
        arb_init(t);
        arb_mul_2exp_si(t, x, 1);
        if (!neg)
            arb_neg(t, t);
        arb_exp(t, t, prec + 4);
        arb_add_ui(y, t, 1, prec + 4);
        arb_div(y, t, y, prec + 4);
        arb_mul_2exp_si(y, y, 1);
        arb_sub_ui(y, y, 1, prec);
        if (!neg)
            arb_neg(y, y);
        if (mag_cmp_2exp_si(arb_radref(y), 0) > 0)
            arb_zero_pm_one(y);
        arb_clear(t);
        return;
    }

    {
        /* the precision bounded by the accuracy of the input, plus the
           bits the derivative sech^2 m < 4 e^(-2 |m|) gains for |m| >= 1:
           2 |m| log2(e) - 2 (|m| in double arithmetic, rounded down) */
        slong e = ARF_EXP(arb_midref(x)), gain = 0;
        if (COEFF_IS_MPZ(e) || e < 30)
        {
            if (!COEFF_IS_MPZ(e) && e >= 1)
            {
                double dm = arf_get_d(arb_midref(x), ARF_RND_DOWN);
                gain = FLINT_MAX((slong) (2.8853900817779268 * ((dm < 0) ? -dm : dm)) - 2, 0);
            }
            prec = FLINT_MAX(FLINT_MIN(prec, 20 - radexp + gain), 2);
        }
    }


    if (radexp >= -9)
    {
        /* the endpoints (rounded outward at prec + 10 bits), the hull
           of the values (evaluated with a few guard bits at the reduced
           precision: a sharp result costs one rounding) */
        arf_t a, b;
        arb_t u;
        slong maglim = FLINT_MAX(128, 2 * prec);
        arf_init(a);
        arf_init(b);
        arb_init(u);
        arb_get_interval_arf(a, b, x, prec + 10);
        _arb_sinh_cosh_arf(NULL, NULL, u, NULL, a, prec + 4, maglim);
        _arb_sinh_cosh_arf(NULL, NULL, y, NULL, b, prec + 4, maglim);
        /* the hull [lower end of tanh(a), upper end of tanh(b)] */
        arf_set_mag(a, arb_radref(u));
        arf_sub(a, arb_midref(u), a, prec + 4, ARF_RND_FLOOR);
        arf_set_mag(b, arb_radref(y));
        arf_add(b, arb_midref(y), b, prec + 4, ARF_RND_CEIL);
        arb_set_interval_arf(y, a, b, prec + 4);
        arb_set_round(y, y, prec);
        if (mag_cmp_2exp_si(arb_radref(y), 0) > 0)
            arb_zero_pm_one(y);
        arf_clear(a);
        arf_clear(b);
        arb_clear(u);
        return;
    }

    mag_init_set(r, arb_radref(x));
    mag_init(t);
    mag_init(u);
    _arb_sinh_cosh_arf(NULL, NULL, y, NULL, arb_midref(x), prec, FLINT_MAX(128, 2 * prec));

    /* the derivative bound: 1 - |tanh m| from above (in arf, as a mag
       lower bound of |tanh m| near 1 would lose it to MAG_BITS) */
    arb_get_mag_lower(t, y);
    if (mag_cmp_2exp_si(r, -1) <= 0 && mag_cmp_2exp_si(t, -1) >= 0)
    {
        arf_t d;
        arf_init(d);
        arf_abs(d, arb_midref(y));
        arf_sub_ui(d, d, 1, MAG_BITS, ARF_RND_DOWN);    /* <= |tanh m| - 1 */
        arf_neg(d, d);
        arf_get_mag(t, d);
        mag_add(t, t, arb_radref(y));
        arf_clear(d);
        mag_mul_2exp_si(u, r, 2);
        mag_add_ui(u, u, 1);
        mag_mul(t, t, u);
        mag_mul_2exp_si(t, t, 1);
        if (mag_cmp_2exp_si(t, 0) > 0)
            mag_one(t);
    }
    else
        mag_one(t);
    mag_mul(t, t, r);

    mag_add(arb_radref(y), arb_radref(y), t);
    if (mag_cmp_2exp_si(arb_radref(y), 0) > 0)
        arb_zero_pm_one(y);
    mag_clear(r);
    mag_clear(t);
    mag_clear(u);
}
