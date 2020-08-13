/*
    Copyright (C) 2018 William Hart

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
_fmpz_poly_sqrt_series(fmpz * res, const fmpz * poly, slong len, slong n)
{
    int result;
    slong i;
    fmpz * rev;
    
    while (len > 0 && n > 0 && fmpz_is_zero(poly))
    {
        if (len > 1 && !fmpz_is_zero(poly + 1))
            return 0;

        fmpz_zero(res);
        fmpz_zero(res + n - 1);

        poly += 2;
        len -= 2;
        n -= 2;
        res++;
    }

    if (len <= 0)
    {
        for (i = 0; i < n; i++)
            fmpz_zero(res + i);
        
        return 1;
    }

    if (n <= 0)
       return 1;

    rev = _fmpz_vec_init(2*n - 1);

    _fmpz_poly_reverse(rev, poly, FLINT_MIN(2*n - 1, len), 2*n - 1);
    result = _fmpz_poly_sqrt_divconquer(res, rev, 2*n - 1, 0);

    if (result)
        _fmpz_poly_reverse(res, res, n, n);

    _fmpz_vec_clear(rev, 2*n - 1);

    return result;
}

int
fmpz_poly_sqrt_series(fmpz_poly_t b, const fmpz_poly_t a, slong n)
{
    slong len = a->length;
    int result;

    if (n == 0 || len == 0)
    {
        fmpz_poly_zero(b);
        return 1;
    }

    if (b == a)
    {
        fmpz_poly_t tmp;
        fmpz_poly_init(tmp);
        result = fmpz_poly_sqrt_series(tmp, a, n);
        fmpz_poly_swap(b, tmp);
        fmpz_poly_clear(tmp);
        return result;
    }

    fmpz_poly_fit_length(b, n);
    _fmpz_poly_set_length(b, n);

    result = _fmpz_poly_sqrt_series(b->coeffs, a->coeffs, len, n);
    if (result)
        _fmpz_poly_normalise(b);
    else
        _fmpz_poly_set_length(b, 0);

    return result;
}
