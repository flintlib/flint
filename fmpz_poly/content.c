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
_fmpz_poly_content(fmpz_t res, const fmpz * poly, slong len)
{
    fmpz_zero(res);
    while (len--)
        fmpz_gcd(res, res, poly + len);
}

void
fmpz_poly_content(fmpz_t res, const fmpz_poly_t poly)
{
    fmpz_t t;
    fmpz_init(t);
    _fmpz_poly_content(t, poly->coeffs, poly->length);
    fmpz_swap(res, t);
    fmpz_clear(t);
}
