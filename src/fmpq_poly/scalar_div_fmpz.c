/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void _fmpq_poly_scalar_div_fmpz(fmpz * rpoly, fmpz_t rden, const fmpz * poly,
                                const fmpz_t den, slong len, const fmpz_t c)
{
    if (fmpz_is_one(c))
    {
        if (rpoly != poly)
        {
            _fmpz_vec_set(rpoly, poly, len);
            fmpz_set(rden, den);
        }
    }
    else if (*c == WORD(-1))
    {
        _fmpz_vec_neg(rpoly, poly, len);
        fmpz_set(rden, den);
    }
    else
    {
        fmpz_t d;
        fmpz_init(d);
        _fmpz_vec_content_chained(d, poly, len, c);

        if (fmpz_sgn(c) < 0)
            fmpz_neg(d, d);
        _fmpz_vec_scalar_divexact_fmpz(rpoly, poly, len, d);
        fmpz_divexact(d, c, d);
        fmpz_mul(rden, den, d);

        fmpz_clear(d);
    }
}

void fmpq_poly_scalar_div_fmpz(fmpq_poly_t rop, const fmpq_poly_t op, const fmpz_t c)
{
    if (*c == WORD(0))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_poly_scalar_div_fmpz). Division by zero.\n");
    }

    if (fmpq_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
        return;
    }

    fmpq_poly_fit_length(rop, op->length);
    _fmpq_poly_set_length(rop, op->length);

    _fmpq_poly_scalar_div_fmpz(rop->coeffs, rop->den,
                               op->coeffs, op->den, op->length, c);
}

