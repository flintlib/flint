/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "padic_poly.h"

void padic_poly_one(padic_poly_t poly)
{
    if (padic_poly_prec(poly) > 0)
    {
        padic_poly_fit_length(poly, 1);
        fmpz_one(poly->coeffs);
        _padic_poly_set_length(poly, 1);
        poly->val = 0;
    }
    else
    {
        padic_poly_zero(poly);
    }
}
