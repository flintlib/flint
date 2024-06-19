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
#include "fmpq_poly.h"

void
fmpq_poly_get_numerator(fmpz_poly_t res, const fmpq_poly_t poly)
{
    fmpz_poly_fit_length(res, poly->length);
    _fmpz_vec_set(res->coeffs, poly->coeffs, poly->length);
    _fmpz_poly_set_length(res, poly->length);
}

void
fmpq_poly_get_denominator(fmpz_t den, const fmpq_poly_t poly)
{
   fmpz_set(den, fmpq_poly_denref(poly));
}
