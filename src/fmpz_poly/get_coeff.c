/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_poly.h"

void
fmpz_poly_get_coeff_fmpz(fmpz_t x, const fmpz_poly_t poly, slong n)
{
    if (n < poly->length)
        fmpz_set(x, poly->coeffs + n);
    else
        fmpz_zero(x);
}

slong
fmpz_poly_get_coeff_si(const fmpz_poly_t poly, slong n)
{
    return (n < poly->length) ? fmpz_get_si(poly->coeffs + n) : WORD(0);
}

ulong
fmpz_poly_get_coeff_ui(const fmpz_poly_t poly, slong n)
{
    return (n < poly->length) ? fmpz_get_ui(poly->coeffs + n) : UWORD(0);
}
