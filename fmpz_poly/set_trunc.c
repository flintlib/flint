/*
    Copyright (C) 2014 Fredrik Johansson

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
fmpz_poly_set_trunc(fmpz_poly_t res, const fmpz_poly_t poly, slong n)
{
    if (poly == res)
    {
        fmpz_poly_truncate(res, n);
    }
    else
    {
        slong rlen;

        rlen = FLINT_MIN(n, poly->length);
        while (rlen > 0 && fmpz_is_zero(poly->coeffs + rlen - 1))
            rlen--;

        fmpz_poly_fit_length(res, rlen);
        _fmpz_vec_set(res->coeffs, poly->coeffs, rlen);
        _fmpz_poly_set_length(res, rlen);
    }
}

