/*
    Copyright (C) 2011 Sebastian Pancratz

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
#include "fmpq_poly.h"

void fmpq_poly_reverse(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
{
    slong len = FLINT_MIN(n, poly->length);

    if (len == 0)
    {
        fmpq_poly_zero(res);
        return;
    }

    fmpq_poly_fit_length(res, n);
    _fmpz_poly_reverse(res->coeffs, poly->coeffs, len, n);
    fmpz_set(res->den, poly->den);
    _fmpq_poly_set_length(res, n);

    fmpq_poly_canonicalise(res);
}

