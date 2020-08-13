/*
    Copyright (C) 2016 William Hart

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
fmpz_poly_scalar_abs(fmpz_poly_t res, const fmpz_poly_t poly)
{
    fmpz_poly_fit_length(res, poly->length);

    _fmpz_vec_scalar_abs(res->coeffs, poly->coeffs, poly->length);

    _fmpz_poly_set_length(res, poly->length);
}
