/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

static void
acb_rising_get_mag2_right(mag_t bound, const arb_t a, const arb_t b, ulong n)
{
    mag_t t, u;
    ulong k;

    mag_init(t);
    mag_init(u);

    arb_get_mag(t, a);
    arb_get_mag(u, b);

    mag_mul(bound, t, t);
    mag_addmul(bound, u, u);
    mag_set(u, bound);
    mag_mul_2exp_si(t, t, 1);

    for (k = 1; k < n; k++)
    {
        mag_add_ui(u, u, 2 * k - 1);
        mag_add(u, u, t);
        mag_mul(bound, bound, u);
    }

    mag_clear(t);
    mag_clear(u);
}

void
acb_rising_ui_get_mag(mag_t bound, const acb_t s, ulong n)
{
    if (n == 0)
    {
        mag_one(bound);
        return;
    }

    if (n == 1)
    {
        acb_get_mag(bound, s);
        return;
    }

    if (!acb_is_finite(s))
    {
        mag_inf(bound);
        return;
    }

    if (arf_sgn(arb_midref(acb_realref(s))) >= 0)
    {
        acb_rising_get_mag2_right(bound, acb_realref(s), acb_imagref(s), n);
    }
    else
    {
        arb_t a;
        slong k;
        mag_t bound2, t, u;

        arb_init(a);
        mag_init(bound2);
        mag_init(t);
        mag_init(u);

        arb_get_mag(u, acb_imagref(s));
        mag_mul(u, u, u);
        mag_one(bound);

        for (k = 0; k < n; k++)
        {
            arb_add_ui(a, acb_realref(s), k, MAG_BITS);

            if (arf_sgn(arb_midref(a)) >= 0)
            {
                acb_rising_get_mag2_right(bound2, a, acb_imagref(s), n - k);
                mag_mul(bound, bound, bound2);
                break;
            }
            else
            {
                arb_get_mag(t, a);
                mag_mul(t, t, t);
                mag_add(t, t, u);
                mag_mul(bound, bound, t);
            }
        }

        arb_clear(a);
        mag_clear(bound2);
        mag_clear(t);
        mag_clear(u);
    }

    mag_sqrt(bound, bound);
}

