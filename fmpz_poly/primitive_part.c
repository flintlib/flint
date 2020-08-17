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
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_primitive_part(fmpz * res, const fmpz * poly, slong len)
{
    fmpz_t x;
    fmpz_init(x);
    _fmpz_poly_content(x, poly, len);
    if (fmpz_sgn(poly + (len - 1)) < 0)
        fmpz_neg(x, x);
    _fmpz_vec_scalar_divexact_fmpz(res, poly, len, x);
    fmpz_clear(x);
}

void
fmpz_poly_primitive_part(fmpz_poly_t res, const fmpz_poly_t poly)
{
    slong len = poly->length;
    if (len == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    fmpz_poly_fit_length(res, len);
    _fmpz_poly_set_length(res, len);
    _fmpz_poly_primitive_part(res->coeffs, poly->coeffs, len);
}
