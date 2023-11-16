/*
    Copyright (C) 2015 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
_arb_sinc_derivative_bound(mag_t d, const arb_t x)
{
    /* |f'(x)| < min(arb_get_mag(x), 1) / 2 */
    mag_t r, one;
    mag_init(r);
    mag_init(one);
    arb_get_mag(r, x);
    mag_one(one);
    mag_min(d, r, one);
    mag_mul_2exp_si(d, d, -1);
    mag_clear(r);
    mag_clear(one);
}

void
_arb_sinc_direct(arb_t z, const arb_t x, slong prec)
{
    /* z = sin(x) / x */
    slong wp;
    arb_t y;
    wp = prec + 2;
    arb_init(y);
    arb_sin(y, x, wp);
    arb_div(z, y, x, prec);
    arb_clear(y);
}

void
arb_sinc(arb_t z, const arb_t x, slong prec)
{
    mag_t c, r;
    mag_init(c);
    mag_init(r);
    mag_set_ui_2exp_si(c, 5, -1);
    arb_get_mag_lower(r, x);
    if (mag_cmp(c, r) < 0)
    {
        /* x is not near the origin */
        _arb_sinc_direct(z, x, prec);
    }
    else if (mag_cmp_2exp_si(arb_radref(x), 1) < 0)
    {
        /* determine error magnitude using the derivative bound */
        if (arb_is_exact(x))
        {
            mag_zero(c);
        }
        else
        {
            _arb_sinc_derivative_bound(r, x);
            mag_mul(c, arb_radref(x), r);
        }

        /* evaluate sinc at the midpoint of x */
        if (arf_is_zero(arb_midref(x)))
        {
            arb_one(z);
        }
        else
        {
            arb_get_mid_arb(z, x);
            _arb_sinc_direct(z, z, prec);
        }

        /* add the error */
        mag_add(arb_radref(z), arb_radref(z), c);
    }
    else
    {
        /* x has a large radius and includes points near the origin */
        arf_zero(arb_midref(z));
        mag_one(arb_radref(z));
    }

    mag_clear(c);
    mag_clear(r);
}
