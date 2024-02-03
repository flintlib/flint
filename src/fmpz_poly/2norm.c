/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_2norm(fmpz_t res, const fmpz * poly, slong len)
{
    _fmpz_vec_dot(res, poly, poly, len);
    fmpz_sqrt(res, res);
}

void
fmpz_poly_2norm(fmpz_t res, const fmpz_poly_t poly)
{
    _fmpz_poly_2norm(res, poly->coeffs, poly->length);
}
