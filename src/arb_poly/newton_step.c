/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

int
_arb_poly_newton_step(arb_t xnew, arb_srcptr poly, slong len,
    const arb_t x,
    const arb_t convergence_interval,
    const arf_t convergence_factor, slong prec)
{
    arf_t err;
    arb_t t, u, v;
    int result;

    arf_init(err);
    arb_init(t);
    arb_init(u);
    arb_init(v);

    arf_set_mag(err, arb_radref(x));
    arf_mul(err, err, err, MAG_BITS, ARF_RND_UP);
    arf_mul(err, err, convergence_factor, MAG_BITS, ARF_RND_UP);

    arf_set(arb_midref(t), arb_midref(x));
    mag_zero(arb_radref(t));

    _arb_poly_evaluate2(u, v, poly, len, t, prec);

    arb_div(u, u, v, prec);
    arb_sub(u, t, u, prec);

    arb_add_error_arf(u, err);

    if (arb_contains(convergence_interval, u) &&
        (mag_cmp(arb_radref(u), arb_radref(x)) < 0))
    {
        arb_swap(xnew, u);
        result = 1;
    }
    else
    {
        arb_set(xnew, x);
        result = 0;
    }

    arb_clear(t);
    arb_clear(u);
    arb_clear(v);
    arf_clear(err);

    return result;
}

