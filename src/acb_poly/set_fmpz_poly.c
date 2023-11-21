/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "acb_poly.h"

void
acb_poly_set_fmpz_poly(acb_poly_t poly, const fmpz_poly_t src, slong prec)
{
    slong i, len = fmpz_poly_length(src);

    acb_poly_fit_length(poly, len);
    _acb_poly_set_length(poly, len);

    for (i = 0; i < len; i++)
        acb_set_round_fmpz(poly->coeffs + i, src->coeffs + i, prec);
}
