/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_calc.h"

int arb_calc_newton_step(arb_t xnew, arb_calc_func_t func,
    void * param, const arb_t x, const arb_t conv_region,
    const arf_t conv_factor, slong prec)
{
    mag_t err, v;
    arb_t t;
    arb_struct u[2];
    int result;

    mag_init(err);
    mag_init(v);
    arb_init(t);
    arb_init(u + 0);
    arb_init(u + 1);

    mag_mul(err, arb_radref(x), arb_radref(x));
    arf_get_mag(v, conv_factor);
    mag_mul(err, err, v);

    arf_set(arb_midref(t), arb_midref(x));
    mag_zero(arb_radref(t));

    func(u, t, param, 2, prec);

    arb_div(u, u, u + 1, prec);
    arb_sub(u, t, u, prec);

    mag_add(arb_radref(u), arb_radref(u), err);

    if (arb_contains(conv_region, u) &&
        (mag_cmp(arb_radref(u), arb_radref(x)) < 0))
    {
        arb_swap(xnew, u);
        result = ARB_CALC_SUCCESS;
    }
    else
    {
        arb_set(xnew, x);
        result = ARB_CALC_NO_CONVERGENCE;
    }

    arb_clear(t);
    arb_clear(u);
    arb_clear(u + 1);
    mag_clear(err);
    mag_clear(v);

    return result;
}

