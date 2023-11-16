/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

static void
arb_sqrt1pm1_tiny(arb_t r, const arb_t z, slong prec)
{
    mag_t b, c;
    arb_t t;

    mag_init(b);
    mag_init(c);
    arb_init(t);

    /* if |z| < 1, then |(sqrt(1+z)-1) - (z/2-z^2/8)| <= |z|^3/(1-|z|)/16 */
    arb_get_mag(b, z);
    mag_one(c);
    mag_sub_lower(c, c, b);
    mag_pow_ui(b, b, 3);
    mag_div(b, b, c);
    mag_mul_2exp_si(b, b, -4);

    arb_mul(t, z, z, prec);
    arb_mul_2exp_si(t, t, -2);
    arb_sub(r, z, t, prec);
    arb_mul_2exp_si(r, r, -1);

    if (mag_is_finite(b))
        arb_add_error_mag(r, b);
    else
        arb_indeterminate(r);

    mag_clear(b);
    mag_clear(c);
    arb_clear(t);
}

void
arb_sqrt1pm1(arb_t r, const arb_t z, slong prec)
{
    slong magz, wp;

    if (arb_is_zero(z))
    {
        arb_zero(r);
        return;
    }

    magz = arf_abs_bound_lt_2exp_si(arb_midref(z));

    if (magz < -prec)
    {
        arb_sqrt1pm1_tiny(r, z, prec);
    }
    else
    {
        if (magz < 0)
            wp = prec + (-magz) + 4;
        else
            wp = prec + 4;

        arb_add_ui(r, z, 1, wp);
        arb_sqrt(r, r, wp);
        arb_sub_ui(r, r, 1, wp);
    }
}

