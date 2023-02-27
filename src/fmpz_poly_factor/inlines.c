/*
    Copyright (C) 2017 Tommy Hofmann 

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define FMPZ_POLY_FACTOR_INLINES_C

#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#undef ulong
#include <gmp.h>
#include "flint.h"
#include "fmpz_poly.h"

void fmpz_poly_factor_get_fmpz_poly(fmpz_poly_t z, const fmpz_poly_factor_t F, slong i)
{
    fmpz_poly_set(z, F->p + i);
}

void fmpz_poly_factor_get_fmpz(fmpz_t z, const fmpz_poly_factor_t F)
{
    fmpz_set(z, &(F->c));
}
