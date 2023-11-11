/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "padic_poly.h"

int padic_poly_debug(const padic_poly_t poly)
{
    flint_printf("(alloc = %wd, length = %wd, val = %wd, N = %wd, vec = ", poly->alloc, poly->length, poly->val, poly->N);
    if (poly->coeffs)
    {
        flint_printf("{");
        _fmpz_vec_print(poly->coeffs, poly->alloc);
        flint_printf("}");
    }
    else
    {
        flint_printf("NULL");
    }
    flint_printf(")");

    return 1;
}
