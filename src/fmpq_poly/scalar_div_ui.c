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

void _fmpq_poly_scalar_div_ui(fmpz * rpoly, fmpz_t rden, const fmpz * poly,
                              const fmpz_t den, slong len, ulong c)
{
    if (c == UWORD(1))
    {
        if (rpoly != poly)
            _fmpz_vec_set(rpoly, poly, len);
        fmpz_set(rden, den);
    }
    else
    {
        fmpz_t d, fc;
        ulong ud;
        fmpz_init(d);
        fmpz_init(fc);
        fmpz_set_ui(fc, c);
        _fmpz_vec_content_chained(d, poly, len, fc);

        ud = fmpz_get_ui(d);  /* gcd of d and c fits into a ulong */

        _fmpz_vec_scalar_divexact_ui(rpoly, poly, len, ud);
        fmpz_mul_ui(rden, den, c / ud);

        fmpz_clear(d);
        fmpz_clear(fc);
    }
}

void fmpq_poly_scalar_div_ui(fmpq_poly_t rop, const fmpq_poly_t op, ulong c)
{
    if (c == UWORD(0))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_poly_scalar_div_ui). Division by zero.\n");
    }

    if (fmpq_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
        return;
    }

    fmpq_poly_fit_length(rop, op->length);
    _fmpq_poly_set_length(rop, op->length);

    _fmpq_poly_scalar_div_ui(rop->coeffs, rop->den,
                             op->coeffs, op->den, op->length, c);
}

