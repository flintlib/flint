/*
    Copyright (C) 2010 Sebastian Pancratz
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

void _fmpq_poly_derivative(fmpz * rpoly, fmpz_t rden, 
                           const fmpz * poly, const fmpz_t den, slong len)
{
    _fmpz_poly_derivative(rpoly, poly, len);
    fmpz_set(rden, den);
    _fmpq_poly_canonicalise(rpoly, rden, len - 1);
}

void fmpq_poly_derivative(fmpq_poly_t res, const fmpq_poly_t poly)
{
    slong len = poly->length;
    if (len < 2)
    {
        fmpq_poly_zero(res);
        return;
    }

    fmpq_poly_fit_length(res, len - 1);
    _fmpq_poly_derivative(res->coeffs, res->den, poly->coeffs, poly->den, len);
    _fmpq_poly_set_length(res, len - 1);
}
