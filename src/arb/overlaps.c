/*
    Copyright (C) 2012, 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

/* decide |xm - ym| <= xr + yr */
int
arb_overlaps(const arb_t x, const arb_t y)
{
    arf_t t;
    mag_t a, b;
    int inexact, result;

    if (!arb_is_finite(x) || !arb_is_finite(y))
    {
        /* special cases: positive and negative infinity */
        if (arf_is_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)))
        {
            if (mag_is_finite(arb_radref(y)) && (arf_is_finite(arb_midref(y)) || (arf_is_inf(arb_midref(y)) && !arf_equal(arb_midref(x), arb_midref(y)))))
                return 0;
        }

        if (arf_is_inf(arb_midref(y)) && mag_is_finite(arb_radref(y)))
        {
            if (mag_is_finite(arb_radref(x)) && (arf_is_finite(arb_midref(x)) || (arf_is_inf(arb_midref(x)) && !arf_equal(arb_midref(x), arb_midref(y)))))
                return 0;
        }

        return 1;
    }

    if (arf_equal(arb_midref(x), arb_midref(y)))
        return 1;

    arf_init(t);
    mag_init(a);
    mag_init(b);

    /* t = lower bound for |xm - ym| */
    inexact = arf_sub(t, arb_midref(x), arb_midref(y), 2 * MAG_BITS, ARF_RND_DOWN);

    /* u = lower bound for |xm - ym|; v = upper bound for xr + yr */
    arf_get_mag_lower(a, t);
    mag_add(b, arb_radref(x), arb_radref(y));

    if (mag_cmp(a, b) > 0)
    {
        result = 0;
    }
    else
    {
        /* a = upper bound for |xm - ym|; b = lower bound for xr + yr */
        arf_get_mag(a, t);

        if (inexact)
        {
            MAG_MAN(a)++;
            MAG_ADJUST_ONE_TOO_LARGE(a);
        }

        mag_add_lower(b, arb_radref(x), arb_radref(y));

        if (mag_cmp(a, b) < 0)
        {
            result = 1;
        }
        else
        {
            arf_struct u[4];

            if (arf_cmp(arb_midref(x), arb_midref(y)) >= 0)
            {
                arf_init_set_shallow(u + 0, arb_midref(x));
                arf_init_neg_shallow(u + 1, arb_midref(y));
            }
            else
            {
                arf_init_neg_shallow(u + 0, arb_midref(x));
                arf_init_set_shallow(u + 1, arb_midref(y));
            }

            arf_init_neg_mag_shallow(u + 2, arb_radref(x));
            arf_init_neg_mag_shallow(u + 3, arb_radref(y));

            arf_sum(t, u, 4, MAG_BITS, ARF_RND_DOWN);
            result = arf_sgn(t) <= 0;
        }
    }

    arf_clear(t);
    mag_clear(a);
    mag_clear(b);

    return result;
}
