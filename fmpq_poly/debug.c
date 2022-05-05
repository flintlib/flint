/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "flint-impl.h"
#include "fmpz.h"

int fmpq_poly_debug(const fmpq_poly_t poly)
{
    slong i;

    printf("{alloc: " WORD_FMT "d, length: " WORD_FMT "d, coeffs:", poly->alloc, poly->length);
    for (i = 0; i < poly->alloc; i++)
    {
        printf(" ");
        fmpz_print(poly->coeffs + i);
    }
    printf(", den: ");
    fmpz_print(poly->den);
    printf("}");

    return 1;
}

