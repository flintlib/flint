/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

slong
_acb_poly_validate_roots(acb_ptr roots,
        acb_srcptr poly, slong len, slong prec)
{
    slong i, j, deg;
    slong isolated, nonisolated, total_isolated;
    acb_ptr deriv;
    acb_ptr tmp;
    int *overlap;

    deg = len - 1;

    deriv = _acb_vec_init(deg);
    overlap = flint_calloc(deg, sizeof(int));
    tmp = flint_malloc(sizeof(acb_struct) * deg);

    _acb_poly_derivative(deriv, poly, len, prec);

    /* compute an inclusion interval for each point */
    for (i = 0; i < deg; i++)
    {
        _acb_poly_root_inclusion(roots + i, roots + i,
            poly, deriv, len, prec);
    }

    /* find which points do not overlap with any other points */
    for (i = 0; i < deg; i++)
    {
        for (j = i + 1; j < deg; j++)
        {
            if (acb_overlaps(roots + i, roots + j))
            {
                overlap[i] = overlap[j] = 1;
            }
        }
    }

    /* count and move all isolated roots to the front of the array */
    total_isolated = 0;
    for (i = 0; i < deg; i++)
        total_isolated += (overlap[i] == 0);

    for (i = 0; i < deg; i++)
        tmp[i] = roots[i];

    isolated = 0;
    nonisolated = 0;
    for (i = 0; i < deg; i++)
    {
        if (overlap[i] == 0)
        {
            roots[isolated] = tmp[i];
            isolated++;
        }
        else
        {
            roots[total_isolated + nonisolated] = tmp[i];
            nonisolated++;
        }
    }

    _acb_vec_clear(deriv, deg);
    flint_free(tmp);
    flint_free(overlap);

    return isolated;
}

