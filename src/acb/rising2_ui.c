/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

void
acb_rising2_ui(acb_t u, acb_t v, const acb_t x, ulong n, slong prec)
{
    if (x == u || x == v)
    {
        acb_t t;
        acb_init(t);
        acb_set(t, x);
        acb_rising2_ui(u, v, t, n, prec);
        acb_clear(t);
    }
    else
    {
        acb_struct tmp[2];

        tmp[0] = *u;
        tmp[1] = *v;

        acb_hypgeom_rising_ui_jet(tmp, x, n, 2, prec);

        *u = tmp[0];
        *v = tmp[1];
    }
}

