/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb_hypgeom.h"

void
arb_rising2_ui(arb_t u, arb_t v, const arb_t x, ulong n, slong prec)
{
    if (x == u || x == v)
    {
        arb_t t;
        arb_init(t);
        arb_set(t, x);
        arb_rising2_ui(u, v, t, n, prec);
        arb_clear(t);
    }
    else
    {
        arb_struct tmp[2];

        tmp[0] = *u;
        tmp[1] = *v;

        arb_hypgeom_rising_ui_jet(tmp, x, n, 2, prec);

        *u = tmp[0];
        *v = tmp[1];
    }
}

