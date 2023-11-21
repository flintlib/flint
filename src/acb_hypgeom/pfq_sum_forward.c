/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_pfq_sum_forward(acb_t s, acb_t t,
    acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, slong n, slong prec)
{
    acb_t u, v;
    slong k, i;

    acb_init(u);
    acb_init(v);

    acb_zero(s);
    acb_one(t);

    for (k = 0; k < n && !acb_is_zero(t); k++)
    {
        acb_add(s, s, t, prec);

        if (p > 0)
        {
            acb_add_ui(u, a, k, prec);

            for (i = 1; i < p; i++)
            {
                acb_add_ui(v, a + i, k, prec);
                acb_mul(u, u, v, prec);
            }

            acb_mul(t, t, u, prec);
        }

        if (q > 0)
        {
            acb_add_ui(u, b, k, prec);

            for (i = 1; i < q; i++)
            {
                acb_add_ui(v, b + i, k, prec);
                acb_mul(u, u, v, prec);
            }

            acb_div(t, t, u, prec);
        }

        acb_mul(t, t, z, prec);
    }

    acb_clear(u);
    acb_clear(v);
}

