/*
    Copyright (C) 2012, 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

/* decide |t| <= xr; t = xm - y */
int
arb_contains_arf(const arb_t x, const arf_t y)
{
    if (!arb_is_finite(x))
    {
        if (arf_is_nan(arb_midref(x)))
            return 1;

        if (arf_is_nan(y))
            return 0;

        if (arf_is_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)))
        {
            return arf_equal(arb_midref(x), y);
        }

        return 1;
    }
    else if (!arf_is_finite(y))
    {
        return 0;
    }
    else if (arb_is_exact(x))
    {
        return arf_equal(arb_midref(x), y);
    }
    else
    {
        arf_t t;
        arf_struct tmp[3];
        int result, inexact;

        arf_init(t);

        inexact = arf_sub(t, arb_midref(x), y, 2 * MAG_BITS, ARF_RND_DOWN);

        if (!inexact)
        {
            result = arf_cmpabs_mag(t, arb_radref(x)) <= 0;
        }
        else
        {
            mag_t a;
            mag_init(a);

            arf_get_mag_lower(a, t);
            if (mag_cmp(a, arb_radref(x)) > 0)
            {
                result = 0;
            }
            else
            {
                arf_get_mag(a, t);
                if (mag_cmp(a, arb_radref(x)) < 0)
                {
                    result = 1;
                }
                else
                {
                    /* y >= xm - xr  <=>  0 >= xm - xr - y */
                    arf_init_set_shallow(tmp + 0, arb_midref(x));
                    arf_init_neg_mag_shallow(tmp + 1,  arb_radref(x));
                    arf_init_neg_shallow(tmp + 2, y);

                    arf_sum(t, tmp, 3, MAG_BITS, ARF_RND_DOWN);
                    result = (arf_sgn(t) <= 0);

                    if (result)
                    {
                        /* y <= xm + xr  <=>  0 <= xm + xr - y */
                        arf_init_set_mag_shallow(tmp + 1,  arb_radref(x));
                        arf_sum(t, tmp, 3, MAG_BITS, ARF_RND_DOWN);
                        result = (arf_sgn(t) >= 0);
                    }
                }
            }

            mag_clear(a);
        }

        arf_clear(t);

        return result;
    }
}
