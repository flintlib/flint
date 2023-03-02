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
#include "fmpz_poly.h"

void
_fmpz_poly_reverse(fmpz * res, const fmpz * poly, slong len, slong n)
{
    if (res == poly)
    {
        slong i;

        for (i = 0; i < n / 2; i++)
        {
            fmpz t = res[i];
            res[i] = res[n - 1 - i];
            res[n - 1 - i] = t;
        }

        for (i = 0; i < n - len; i++)
            fmpz_zero(res + i);
    }
    else
    {
        slong i;

        for (i = 0; i < n - len; i++)
            fmpz_zero(res + i);

        for (i = 0; i < len; i++)
            fmpz_set(res + (n - len) + i, poly + (len - 1) - i);
    }
}

void
fmpz_poly_reverse(fmpz_poly_t res, const fmpz_poly_t poly, slong n)
{
    slong len = FLINT_MIN(n, poly->length);
    if (len == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    fmpz_poly_fit_length(res, n);

    _fmpz_poly_reverse(res->coeffs, poly->coeffs, len, n);

    _fmpz_poly_set_length(res, n);
    _fmpz_poly_normalise(res);
}
