/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_doublefac_ui(arb_t res, ulong n, slong prec)
{
    if (n % 2 == 0)
    {
        arb_fac_ui(res, n / 2, prec);
        arb_mul_2exp_si(res, res, n / 2);
    }
    else
    {
        arb_t t;
        arb_init(t);
        arb_doublefac_ui(t, n - 1, prec + 5);
        arb_fac_ui(res, n, prec + 5);
        arb_div(res, res, t, prec);
        arb_clear(t);
    }
}

