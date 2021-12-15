/*
    Copyright (C) 2011, 2010 Sebastian Pancratz
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

void _fmpz_poly_derivative(fmpz * rpoly, const fmpz * poly, slong len)
{
    slong i;

    for (i = 1; i < len; i++)
        fmpz_mul_ui(rpoly + (i - 1), poly + i, i);
}

void fmpz_poly_derivative(fmpz_poly_t res, const fmpz_poly_t poly)
{
    const slong len = poly->length;

    if (len < 2)
    {
        fmpz_poly_zero(res);
    }
    else
    {
        fmpz_poly_fit_length(res, len - 1);
        _fmpz_poly_derivative(res->coeffs, poly->coeffs, len);
        _fmpz_poly_set_length(res, len - 1);
    }
}
