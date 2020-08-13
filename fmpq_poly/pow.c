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
#include "fmpq_poly.h"

void _fmpq_poly_pow(fmpz * rpoly, fmpz_t rden, 
                    const fmpz * poly, const fmpz_t den, slong len, ulong e)
{
    _fmpz_poly_pow(rpoly, poly, len, e);
    fmpz_pow_ui(rden, den, e);
}

void fmpq_poly_pow(fmpq_poly_t res, const fmpq_poly_t poly, ulong e)
{
    slong len = poly->length, rlen;

    if (e == 0)
    {
        fmpq_poly_set_ui(res, 1);
        return;
    }
    if (len == 0)
    {
        fmpq_poly_zero(res);
        return;
    }

    rlen = (slong) e * (len - WORD(1)) + WORD(1);

    if (res != poly)
    {
        fmpq_poly_fit_length(res, rlen);
        _fmpq_poly_pow(res->coeffs, res->den, poly->coeffs, poly->den, len, e);
        _fmpq_poly_set_length(res, rlen);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, rlen);
        _fmpq_poly_pow(t->coeffs, t->den, poly->coeffs, poly->den, len, e);
        _fmpq_poly_set_length(t, rlen);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
    }
}

