/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "arb_poly.h"

int
arb_poly_get_unique_fmpz_poly(fmpz_poly_t res, const arb_poly_t src)
{
    int success;

    fmpz_poly_fit_length(res, src->length);
    success = _arb_vec_get_unique_fmpz_vec(res->coeffs, src->coeffs, src->length);
    _fmpz_poly_set_length(res, src->length);
    _fmpz_poly_normalise(res);
    return success;
}

