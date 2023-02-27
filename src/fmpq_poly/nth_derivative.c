/*
    Copyright (C) 2021 Mathieu Gouttenoire

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

void _fmpq_poly_nth_derivative(fmpz * rpoly, fmpz_t rden, 
    const fmpz * poly, const fmpz_t den, ulong n, slong len)
{
    _fmpz_poly_nth_derivative(rpoly, poly, n, len);
    fmpz_set(rden, den);
    _fmpq_poly_canonicalise(rpoly, rden, len - n);
}

void fmpq_poly_nth_derivative(fmpq_poly_t res, const fmpq_poly_t poly, ulong n)
{
    slong len = poly->length;
    if (len <= n)
    {
        fmpq_poly_zero(res);
        return;
    }
    
    fmpq_poly_fit_length(res, len - n);
    if (n == 0)
    {
        fmpq_poly_set(res, poly);
    }
    else if (n == 1)
    {
        _fmpq_poly_derivative(res->coeffs, res->den, poly->coeffs, poly->den, len);
    }
    else
    {
        _fmpq_poly_nth_derivative(res->coeffs, res->den, poly->coeffs, poly->den, n, len);
    }
    _fmpq_poly_set_length(res, len - n);
}
