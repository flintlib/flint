/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void _fmpq_poly_scalar_mul_si(fmpz * rpoly, fmpz_t rden, 
                              const fmpz * poly, const fmpz_t den, slong len, 
                              slong c)
{
    fmpz_t gcd;  /* GCD( den, c ) */

    if (c == 0)
    {
        _fmpz_vec_zero(rpoly, len);
        fmpz_one(rden);
        return;
    }

    fmpz_init(gcd);
    fmpz_set_si(gcd, c);
    fmpz_gcd(gcd, gcd, den);
    if (fmpz_is_one(gcd))
    {
        _fmpz_vec_scalar_mul_si(rpoly, poly, len, c);
        fmpz_set(rden, den);
    }
    else
    {
        if (c > WORD_MIN || fmpz_cmp_ui(gcd, - (ulong) WORD_MIN))
        {
            slong g = fmpz_get_si(gcd);

            _fmpz_vec_scalar_mul_si(rpoly, poly, len, c / g);
            fmpz_divexact_si(rden, den, g);
        }
        else
        {
            _fmpz_vec_neg(rpoly, poly, len);
            fmpz_divexact_ui(rden, den, - (ulong) WORD_MIN);
        }
    }
    fmpz_clear(gcd);
}

void fmpq_poly_scalar_mul_si(fmpq_poly_t rop, const fmpq_poly_t op, slong c)
{
    if (c == 0 || fmpq_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
        return;
    }
    
    fmpq_poly_fit_length(rop, op->length);
    _fmpq_poly_set_length(rop, op->length);
    
    _fmpq_poly_scalar_mul_si(rop->coeffs, rop->den, 
                             op->coeffs, op->den, op->length, c);
}

