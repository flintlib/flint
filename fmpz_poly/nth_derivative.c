/*
    Copyright (C) 2021, Mathieu Gouttenoire

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

void _fmpz_poly_nth_derivative(fmpz * rpoly, const fmpz * poly, ulong n, slong len)
{
    slong i;
    fmpz_t c;
    fmpz_init(c);
    
    fmpz_fac_ui(c, n);
    fmpz_mul(rpoly, poly + n, c);
    for (i = n + 1; i < len; i ++)
    {
        fmpz_divexact_ui(c, c, i - n);
        fmpz_mul_ui(c, c, i);
        fmpz_mul(rpoly + i - n, poly + i, c);
    }
    
    fmpz_clear(c);
}

void fmpz_poly_nth_derivative(fmpz_poly_t res, const fmpz_poly_t poly, ulong n)
{
    const slong len = poly->length;
    
    if (len <= n)
    {
        fmpz_poly_zero(res);
        return;
    }

    fmpz_poly_fit_length(res, len - n);
    if (n == 0)
    {
        fmpz_poly_set(res, poly);
    }
    else if (n == 1)
    {
        _fmpz_poly_derivative(res->coeffs, poly->coeffs, len);
    }
    else
    {
        _fmpz_poly_nth_derivative(res->coeffs, poly->coeffs, n, len);
    }
    _fmpz_poly_set_length(res, len - n);
}
