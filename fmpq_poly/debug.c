/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

int fmpq_poly_debug(const fmpq_poly_t poly)
{
    slong i;

    flint_printf("{alloc: %wd, length: %wd, coeffs:", poly->alloc, poly->length);
    for (i = 0; i < poly->alloc; i++)
    {
        flint_printf(" ");
        fmpz_print(poly->coeffs + i);
    }
    flint_printf(", den: ");
    fmpz_print(poly->den);
    flint_printf("}");

    return 1;
}

