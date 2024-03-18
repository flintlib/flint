/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

ulong fmpz_poly_max_limbs(const fmpz_poly_t poly)
{
    return _fmpz_vec_max_limbs(poly->coeffs, poly->length);
}

slong fmpz_poly_max_bits(const fmpz_poly_t poly)
{
    return _fmpz_vec_max_bits(poly->coeffs, poly->length);
}

void fmpz_poly_height(fmpz_t res, const fmpz_poly_t poly)
{
    _fmpz_vec_height(res, poly->coeffs, poly->length);
}

slong _fmpz_poly_hamming_weight(const fmpz * a, slong len)
{
    slong i, sum = 0;
    for (i = 0; i < len; i++)
        sum += !fmpz_is_zero(a + i);
    return sum;
}
