/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_calc.h"

void arb_calc_newton_conv_factor(arf_t conv_factor,
    arb_calc_func_t func, void * param,
    const arb_t conv_region, slong prec)
{
    arb_struct t[3];

    arb_init(t);
    arb_init(t + 1);
    arb_init(t + 2);

    func(t, conv_region, param, 3, prec);

    arb_div(t, t + 2, t + 1, prec);
    arb_mul_2exp_si(t, t, -1);

    arb_get_abs_ubound_arf(conv_factor, t, prec);

    arb_clear(t);
    arb_clear(t + 1);
    arb_clear(t + 2);
}

