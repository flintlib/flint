/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_2f1_series_direct(acb_poly_t res,
    const acb_poly_t a, const acb_poly_t b,
    const acb_poly_t c, const acb_poly_t z, int regularized, slong len, slong prec)
{
    acb_poly_struct aa[4];

    aa[0] = *a;
    aa[1] = *b;
    aa[2] = *c;
    acb_poly_init(&aa[3]);
    acb_poly_one(&aa[3]);

    acb_hypgeom_pfq_series_direct(res, aa, 2, aa + 2, 2, z,
        regularized, -1, len, prec);

    acb_poly_clear(&aa[3]);
}

