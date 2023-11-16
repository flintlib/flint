/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void
acb_dirichlet_eta(acb_t res, const acb_t s, slong prec)
{
    if (!acb_is_finite(s))
    {
        acb_indeterminate(res);
    }
    else if (arb_contains_si(acb_realref(s), 1) && arb_contains_zero(acb_imagref(s)))
    {
        if (acb_is_one(s))
        {
            arb_const_log2(acb_realref(res), prec);
            arb_zero(acb_imagref(res));
        }
        else
        {
            mag_t m;
            int is_real = acb_is_real(s);
            mag_init(m);

            /* Taylor coefficients at s = 1 are bounded by |c_k| < 1/4^k. */
            acb_sub_ui(res, s, 1, prec);
            acb_get_mag(m, res);
            mag_mul_2exp_si(m, m, -2);
            mag_geom_series(m, m, 1);

            arb_const_log2(acb_realref(res), prec);
            arb_zero(acb_imagref(res));

            if (is_real)
                arb_add_error_mag(acb_realref(res), m);
            else
                acb_add_error_mag(res, m);

            mag_clear(m);
        }
    }
    else
    {
        acb_t t;
        acb_init(t);
        acb_one(t);
        acb_mul_2exp_si(t, t, -1);
        acb_pow(t, t, s, prec);
        acb_mul_2exp_si(t, t, 1);
        acb_sub_ui(t, t, 1, prec);
        acb_neg(t, t);
        acb_zeta(res, s, prec);
        acb_mul(res, res, t, prec);
        acb_clear(t);
    }
}

