/*
    Copyright (C) 2012-2018 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_rsqrt(arb_t z, const arb_t x, slong prec)
{
    mag_t t, u;
    slong acc;
    int inexact;

    if (!arb_is_finite(x) || arf_sgn(arb_midref(x)) <= 0)
    {
        if (arf_is_pos_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)))
            arb_zero(z);
        else
            arb_indeterminate(z);
    }
    else if (mag_is_zero(arb_radref(x)))
    {
        inexact = arf_rsqrt(arb_midref(z), arb_midref(x), prec, ARB_RND);

        if (inexact)
            arf_mag_set_ulp(arb_radref(z), arb_midref(z), prec);
        else
            mag_zero(arb_radref(z));
    }
    else
    {
        mag_init(t);

        arb_get_mag_lower(t, x);

        acc = arb_rel_accuracy_bits(x);
        acc = FLINT_MIN(acc, prec);
        prec = FLINT_MIN(prec, acc + MAG_BITS);
        prec = FLINT_MAX(prec, 2);

        if (acc <= 20)
        {
            if (mag_is_zero(t))
            {
                arb_indeterminate(z);
            }
            else
            {
                mag_init(u);
                arb_get_mag(u, x);

                mag_rsqrt(t, t);
                mag_rsqrt_lower(u, u);

                arb_set_interval_mag(z, u, t, prec);

                mag_clear(u);
            }
        }
        else
        {
            /* error bound: (1/2) (x-r)^(-3/2) * r */
            mag_init(u);

            mag_rsqrt(u, t);
            mag_div(t, u, t);
            mag_mul(t, t, arb_radref(x));
            mag_mul_2exp_si(t, t, -1);

            inexact = arf_rsqrt(arb_midref(z), arb_midref(x), prec, ARB_RND);

            if (inexact)
                arf_mag_add_ulp(arb_radref(z), t, arb_midref(z), prec);
            else
                mag_swap(arb_radref(z), t);

            mag_clear(u);
        }

        mag_clear(t);
    }
}

void
arb_rsqrt_ui(arb_t z, ulong x, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_set_ui(t, x);
    arb_rsqrt(z, t, prec);
    arb_clear(t);
}

