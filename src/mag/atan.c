/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"
#include "arb.h"  /* for atan table */

#ifdef __GNUC__
# define ldexp __builtin_ldexp
#else
# include <math.h>
#endif

static double
mag_atan_d(double x)
{
    double t, u, v;
    int p, q, flip;

    t = x;
    flip = 0;
    if (t > 1.0)
    {
        t = 1.0 / t;
        flip = 1;
    }

    q = ldexp(1.0, ARB_ATAN_TAB1_BITS);
    p = t * (1 << ARB_ATAN_TAB1_BITS);
    if (p == (1 << ARB_ATAN_TAB1_BITS))
        p--;

#if FLINT_BITS == 64
    v = arb_atan_tab1[p][ARB_ATAN_TAB1_LIMBS - 1] * ldexp(1.0, -FLINT_BITS);
#else
    v = arb_atan_tab1[p][ARB_ATAN_TAB1_LIMBS - 1] * ldexp(1.0, -FLINT_BITS)
      + arb_atan_tab1[p][ARB_ATAN_TAB1_LIMBS - 2] * ldexp(1.0, -2 * FLINT_BITS);
#endif

    t = (q * t - p) / (p * t + q);

    /* t - t^3/3 + t^5/5 */
    u = t * t;
    t = (t * (1.0 / 15.0)) * (15.0 - 5.0 * u + 3.0 * (u * u));
    v += t;

    if (flip)
        v = 1.5707963267948966192 - v;

    return v;
}

void
mag_atan(mag_t res, const mag_t x)
{
    if (mag_is_zero(x))
    {
        mag_zero(res);
    }
    else if (mag_cmp_2exp_si(x, MAG_BITS) > 0)
    {
        mag_const_pi(res);
        mag_mul_2exp_si(res, res, -1);
    }
    else if (mag_cmp_2exp_si(x, -(MAG_BITS / 2)) < 0)
    {
        mag_set(res, x);
    }
    else
    {
        double t;

        t = ldexp(MAG_MAN(x), MAG_EXP(x) - MAG_BITS);
        t = mag_atan_d(t);
        t = t * (1.0 + 1e-12);

        mag_set_d(res, t);
    }
}

void
mag_atan_lower(mag_t res, const mag_t x)
{
    if (mag_is_zero(x))
    {
        mag_zero(res);
    }
    else if (mag_is_inf(x))
    {
        mag_const_pi_lower(res);
        mag_mul_2exp_si(res, res, -1);
    }
    else if (mag_cmp_2exp_si(x, MAG_BITS) > 0)
    {
        /* atan(x) > pi/2 - 1/x */
        mag_t t;
        mag_init(t);
        mag_one(t);
        mag_div(t, t, x);
        mag_const_pi_lower(res);
        mag_mul_2exp_si(res, res, -1);
        mag_sub_lower(res, res, t);
        mag_clear(t);
    }
    else if (mag_cmp_2exp_si(x, -MAG_BITS) < 0)
    {
        /* atan(x) > x - x^2 */
        mag_t t;
        mag_init(t);
        mag_mul(t, x, x);
        mag_sub_lower(res, x, t);
        mag_clear(t);
    }
    else
    {
        double t;

        t = ldexp(MAG_MAN(x), MAG_EXP(x) - MAG_BITS);
        t = mag_atan_d(t);
        t = t * (1.0 - 1e-12);

        mag_set_d_lower(res, t);
    }
}

