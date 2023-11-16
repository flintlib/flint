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
arb_asin(arb_t z, const arb_t x, slong prec)
{
    arb_t t;

    if (arb_is_exact(x))
    {
        int side;

        if (arf_is_zero(arb_midref(x)))
        {
            arb_zero(z);
            return;
        }

        side = arf_cmpabs_2exp_si(arb_midref(x), 0);

        /* +/- 1 */
        if (side == 0)
        {
            if (arf_is_one(arb_midref(x)))
            {
                arb_const_pi(z, prec);
            }
            else
            {
                arb_const_pi(z, prec);
                arb_neg(z, z);
            }

            arb_mul_2exp_si(z, z, -1);
            return;
        }
        else if (side > 0)
        {
            arb_indeterminate(z);
            return;
        }
    }

    arb_init(t);

    arb_one(t);
    arb_submul(t, x, x, prec);
    arb_rsqrt(t, t, prec);
    arb_mul(t, x, t, prec);
    arb_atan(z, t, prec);

    arb_clear(t);
}

