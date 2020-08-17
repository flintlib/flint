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
#include "fmpz_poly.h"

void
_fmpz_poly_2norm(fmpz_t res, const fmpz * poly, slong len)
{
    slong i;
    fmpz_zero(res);
    for (i = 0; i < len; i++)
        fmpz_addmul(res, poly + i, poly + i);
    fmpz_sqrt(res, res);
}

void
fmpz_poly_2norm(fmpz_t res, const fmpz_poly_t poly)
{
    _fmpz_poly_2norm(res, poly->coeffs, poly->length);
}
