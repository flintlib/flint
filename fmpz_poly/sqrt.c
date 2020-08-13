/*
    Copyright (C) 2012 Fredrik Johansson

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

int
_fmpz_poly_sqrt(fmpz * res, const fmpz * poly, slong len)
{
    /* KS is heuristic but fast */
    int result = _fmpz_poly_sqrt_KS(res, poly, len);

    /* KS failed, so use fallback */
    if (result == -1)
       result = _fmpz_poly_sqrt_divconquer(res, poly, len, 1);

    return result;
}

int
fmpz_poly_sqrt(fmpz_poly_t b, const fmpz_poly_t a)
{
    slong blen, len = a->length;
    int result;

    if (len % 2 == 0)
    {
        fmpz_poly_zero(b);
        return len == 0;
    }

    if (b == a)
    {
        fmpz_poly_t tmp;
        fmpz_poly_init(tmp);
        result = fmpz_poly_sqrt(tmp, a);
        fmpz_poly_swap(b, tmp);
        fmpz_poly_clear(tmp);
        return result;
    }

    blen = len / 2 + 1;
    fmpz_poly_fit_length(b, blen);
    _fmpz_poly_set_length(b, blen);
    result = _fmpz_poly_sqrt(b->coeffs, a->coeffs, len);
    if (!result)
        _fmpz_poly_set_length(b, 0);
    return result;
}
