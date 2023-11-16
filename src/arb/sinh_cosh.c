/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_sinh_cosh_wide(arb_t s, arb_t c, const arb_t x, slong prec)
{
    mag_t t, u, v, w;

    mag_init(t);
    mag_init(u);
    mag_init(v);
    mag_init(w);

    arb_get_mag_lower(t, x);
    arb_get_mag(u, x);

    if (c != NULL)
    {
        mag_cosh_lower(v, t);
        mag_cosh(w, u);
    }

    if (s != NULL)
    {
        if (mag_is_zero(t))
        {
            arf_get_mag_lower(t, arb_midref(x));
            mag_sub(t, arb_radref(x), t);

            mag_sinh(t, t);
            mag_sinh(u, u);

            if (arf_sgn(arb_midref(x)) > 0)
                arb_set_interval_neg_pos_mag(s, t, u, prec);
            else
                arb_set_interval_neg_pos_mag(s, u, t, prec);
        }
        else
        {
            mag_sinh_lower(t, t);
            mag_sinh(u, u);

            if (arf_sgn(arb_midref(x)) > 0)
            {
                arb_set_interval_mag(s, t, u, prec);
            }
            else
            {
                arb_set_interval_mag(s, t, u, prec);
                arb_neg(s, s);
            }
        }
    }

    if (c != NULL)
        arb_set_interval_mag(c, v, w, prec);

    mag_clear(t);
    mag_clear(u);
    mag_clear(v);
    mag_clear(w);
}

void
arb_cosh(arb_t c, const arb_t x, slong prec)
{
    if (arb_is_zero(x))
    {
        arb_one(c);
    }
    else if (!arb_is_finite(x))
    {
        if (arf_is_nan(arb_midref(x)))
            arb_indeterminate(c);
        else if (arf_is_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)))
            arb_pos_inf(c);
        else
            arb_zero_pm_inf(c);
    }
    else if (mag_cmp_2exp_si(arb_radref(x), -20) > 0 &&
        mag_cmp_2exp_si(arb_radref(x), 10) < 0 &&
        arf_cmpabs_2exp_si(arb_midref(x), 4) < 0)
    {
        arb_sinh_cosh_wide(NULL, c, x, prec);
    }
    else
    {
        arb_t t;
        slong wp = prec + 4;

        /* todo: close to 0, could use derivative to bound
                 propagated error more tightly */
        arb_init(t);
        arb_exp_invexp(c, t, x, wp);
        arb_add(c, c, t, prec);
        arb_mul_2exp_si(c, c, -1);
        arb_clear(t);
    }
}

void
arb_sinh(arb_t s, const arb_t x, slong prec)
{
    if (arb_is_zero(x))
    {
        arb_zero(s);
    }
    else if (!arb_is_finite(x))
    {
        if (arf_is_nan(arb_midref(x)))
            arb_indeterminate(s);
        else if (arf_is_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)))
        {
            arf_set(arb_midref(s), arb_midref(x));
            mag_zero(arb_radref(s));
        }
        else
            arb_zero_pm_inf(s);
    }
    else if (mag_cmp_2exp_si(arb_radref(x), -20) > 0 &&
        mag_cmp_2exp_si(arb_radref(x), 10) < 0 &&
        arf_cmpabs_2exp_si(arb_midref(x), 4) < 0)
    {
        /* for sinh this is not more accurate, but slightly faster */
        arb_sinh_cosh_wide(s, NULL, x, prec);
    }
    else
    {
        arb_t t;
        slong wp = prec + 4;

        arb_init(t);

        if (arf_cmpabs_2exp_si(arb_midref(x), -1) <= 0 &&
               mag_cmp_2exp_si(arb_radref(x), -4) <= 0)
        {
            arb_expm1(s, x, wp);
            arb_add_ui(t, s, 1, wp);
            arb_div(t, s, t, wp);
            arb_add(s, s, t, prec);
        }
        else
        {
            arb_exp_invexp(s, t, x, wp);
            arb_sub(s, s, t, prec);
        }

        arb_mul_2exp_si(s, s, -1);
        arb_clear(t);
    }
}

void
arb_sinh_cosh(arb_t s, arb_t c, const arb_t x, slong prec)
{
    if (arb_is_zero(x))
    {
        arb_zero(s);
        arb_one(c);
    }
    else if (!arb_is_finite(x))
    {
        if (arf_is_nan(arb_midref(x)))
        {
            arb_indeterminate(s);
            arb_indeterminate(c);
        }
        else if (arf_is_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)))
        {
            arf_set(arb_midref(s), arb_midref(x));
            mag_zero(arb_radref(s));
            arf_abs(arb_midref(c), arb_midref(s));
            mag_zero(arb_radref(c));
        }
        else
        {
            arb_zero_pm_inf(s);
            arb_zero_pm_inf(c);
        }
    }
    else if (mag_cmp_2exp_si(arb_radref(x), -20) > 0 &&
        mag_cmp_2exp_si(arb_radref(x), 10) < 0 &&
        arf_cmpabs_2exp_si(arb_midref(x), 4) < 0)
    {
        arb_sinh_cosh_wide(s, c, x, prec);
    }
    else
    {
        arb_t t;
        slong wp = prec + 4;

        arb_init(t);

        if (arf_cmpabs_2exp_si(arb_midref(x), -1) <= 0 &&
            mag_cmp_2exp_si(arb_radref(x), -4) <= 0)
        {
            arb_expm1(s, x, wp);
            arb_add_ui(t, s, 1, wp);
            arb_inv(c, t, wp);
            arb_addmul(s, s, c, prec);
            arb_add(c, c, t, prec);
        }
        else
        {
            arb_exp_invexp(c, t, x, wp);
            arb_sub(s, c, t, prec);
            arb_add(c, c, t, prec);
        }

        arb_mul_2exp_si(s, s, -1);
        arb_mul_2exp_si(c, c, -1);

        arb_clear(t);
    }
}
