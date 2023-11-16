/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_legendre_p_uiui_rec(acb_t res, ulong n, ulong m, const acb_t z, slong prec)
{
    acb_t t, u, v;
    slong k;

    if (!acb_is_finite(z))
    {
        acb_indeterminate(res);
        return;
    }

    if (m > n)
    {
        acb_zero(res);
        return;
    }

    if ((n - m) / 4 > prec)
    {
        acb_indeterminate(res);
        return;
    }

    acb_init(t);
    acb_init(u);
    acb_init(v);

    /* t = p(m,m) = (-1)^m (2m-1)!! */
    if (m == 0)
        arb_one(acb_realref(t));
    else
        arb_doublefac_ui(acb_realref(t), 2 * m - 1, prec);

    if (m % 2)
        arb_neg(acb_realref(t), acb_realref(t));

    if (n > m)
    {
        /* t = p(m+1,m) = z(2m+1)p(m,m), u = p(m,m) */
        acb_mul_ui(u, t, 2 * m + 1, prec);
        acb_mul(u, u, z, prec);
        acb_swap(t, u);

        for (k = m + 2; k <= n; k++)
        {
            /* t, u = ((2*k-1)*z*t - (k+m-1)*u) / (k-m), t */
            acb_mul(v, t, z, prec);
            acb_mul_ui(v, v, 2 * k - 1, prec);
            acb_submul_ui(v, u, k + m - 1, prec);
            acb_div_ui(u, v, k - m, prec);
            acb_swap(t, u);
        }
    }

    acb_set(res, t);

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

