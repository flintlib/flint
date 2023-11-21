/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_lambertw_bound_deriv(mag_t res, const acb_t z, const acb_t ez1, const fmpz_t k)
{
    mag_t t, u, v;

    mag_init(t);
    mag_init(u);
    mag_init(v);

    if (fmpz_is_zero(k))
    {
        acb_get_mag(t, z);

        /* |z| <= 64 */
        if (mag_cmp_2exp_si(t, 6) < 0)
        {
            /* 2.25 / sqrt(|t|) / sqrt(1+|t|), t = |ez+1| */
            acb_get_mag_lower(t, ez1);
            mag_one(u);
            mag_add_lower(u, u, t);
            mag_mul_lower(t, t, u);
            mag_rsqrt(t, t);

            if (arb_is_positive(acb_realref(ez1)))
            {
                mag_mul_ui(t, t, 135);      /* x0.9375 for small improvement */
                mag_mul_2exp_si(t, t, -6);
            }
            else
            {
                mag_mul_ui(t, t, 9);
                mag_mul_2exp_si(t, t, -2);
            }

            mag_set(res, t);
        }
        else
        {
            acb_get_mag_lower(t, z);

            if (mag_cmp_2exp_si(t, 2) >= 0)
            {
                mag_one(u);
                mag_div(res, u, t);
            }
            else   /* unlikely */
            {
                acb_get_mag_lower(u, ez1);
                mag_rsqrt(u, u);
                mag_mul_2exp_si(u, u, -1);
                mag_add_ui(u, u, 1);
                mag_mul_ui(u, u, 3);
                mag_div(res, u, t);
            }
        }
    }
    else if (fmpz_is_pm1(k))
    {
        if (arb_is_nonnegative(acb_realref(z)) ||
            (fmpz_is_one(k) && arb_is_nonnegative(acb_imagref(z))) ||
            (fmpz_equal_si(k, -1) && arb_is_negative(acb_imagref(z))))
        {
            /* (1 + 1/(4+|z|^2))/|z| */
            acb_get_mag_lower(t, z);
            mag_mul_lower(u, t, t);
            mag_set_ui_lower(v, 4);
            mag_add_lower(u, u, v);
            mag_one(v);
            mag_div(u, v, u);
            mag_add(u, u, v);
            mag_div(res, u, t);
        }
        else
        {
            /* (1 + 0.71875/sqrt(|e*z+1|)) / |z| */
            acb_get_mag_lower(t, ez1);
            mag_rsqrt(t, t);
            mag_mul_ui(t, t, 23);
            mag_mul_2exp_si(t, t, -5);
            mag_one(u);
            mag_add(t, t, u);
            acb_get_mag_lower(u, z);
            mag_div(res, t, u);
        }

        mag_clear(t);
        mag_clear(u);
        mag_clear(v);
        return;
    }
    else
    {
        /* |W'(z)|  = |1/z| |W(z)/(1+W(z))|
                   <= |1/z| (2pi) / (2pi-1) */
        mag_set_ui_2exp_si(t, 77, -6);
        acb_get_mag_lower(res, z);
        mag_div(res, t, res);
    }

    mag_clear(t);
    mag_clear(u);
    mag_clear(v);
}

