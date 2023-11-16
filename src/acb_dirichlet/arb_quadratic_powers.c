/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

/* TODO: BSGS can reduce to nv mul */
void
acb_dirichlet_arb_quadratic_powers(arb_ptr v, slong nv, const arb_t x, slong prec)
{
    slong i;
    arb_t dx, x2;
    arb_init(dx);
    arb_init(x2);
    arb_set(dx, x);
    arb_mul(x2, x, x, prec);
    for (i = 0; i < nv; i++)
    {
        if (i == 0)
            arb_one(v + i);
        else if (i == 1)
            arb_set_round(v + i, x, prec);
        else
        {
            arb_mul(dx, dx, x2, prec);
            arb_mul(v + i, v + i - 1, dx, prec);
        }
    }
    arb_clear(dx);
    arb_clear(x2);
}
