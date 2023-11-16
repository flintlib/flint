/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_pfq_series_sum(acb_poly_t s, acb_poly_t t,
    const acb_poly_struct * a, slong p,
    const acb_poly_struct * b, slong q,
    const acb_poly_t z, int regularized,
    slong n, slong len, slong prec)
{
    slong abits, zbits, i, j, cb;

    if (n < 4)
    {
        acb_hypgeom_pfq_series_sum_forward(s, t, a, p, b, q, z,
            regularized, n, len, prec);
        return;
    }

    abits = 0;
    zbits = 0;

    for (i = 0; i < p; i++)
    {
        for (j = 0; j < FLINT_MIN(acb_poly_length(a + i), n); j++)
        {
            cb = acb_bits((a + i)->coeffs + j);
            abits = FLINT_MAX(abits, cb);
        }
    }

    for (i = 0; i < q; i++)
    {
        for (j = 0; j < FLINT_MIN(acb_poly_length(b + i), n); j++)
        {
            cb = acb_bits((b + i)->coeffs + j);
            abits = FLINT_MAX(abits, cb);
        }
    }

    for (j = 0; j < FLINT_MIN(acb_poly_length(z), n); j++)
    {
        cb = acb_bits(z->coeffs + j);
        zbits = FLINT_MAX(zbits, cb);
    }

    /* Prefer RS with small coefficients in parameters, large coefficients
       in z. TODO: tune for larger len? */
    if (len <= 5 && prec > 900 && abits < prec * 0.1 && zbits > prec * 0.1)
    {
        acb_hypgeom_pfq_series_sum_rs(s, t, a, p, b, q, z,
            regularized, n, len, prec);
        return;
    }

    /* Prefer BS with small coefficients and high precision, or when
       computing derivatives of high order. */
    if ((abits < prec * 0.1 && zbits < prec * 0.1 && prec > 600) || len > 20)
    {
        acb_hypgeom_pfq_series_sum_bs(s, t, a, p, b, q, z,
            regularized, n, len, prec);
        return;
    }

    /* TODO: also use bs here when n is large enough, for better
       numerical stability? */

    acb_hypgeom_pfq_series_sum_forward(s, t, a, p, b, q, z,
        regularized, n, len, prec);
}

