/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_poly_mini.h"

void
fmpz_poly_scalar_mod_fmpz(fmpz_poly_t poly1, const fmpz_poly_t poly2, const fmpz_t x)
{
    if (poly2->length == 0)
    {
        fmpz_poly_zero(poly1);
    }
    else
    {
        fmpz_poly_fit_length(poly1, poly2->length);
        _fmpz_vec_scalar_mod_fmpz(poly1->coeffs, poly2->coeffs, poly2->length, x);
        _fmpz_poly_set_length(poly1, poly2->length);
        _fmpz_poly_normalise(poly1);
    }
}
