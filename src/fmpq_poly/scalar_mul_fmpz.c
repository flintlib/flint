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
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void _fmpq_poly_scalar_mul_fmpz(fmpz * rpoly, fmpz_t rden, 
                                const fmpz * poly, const fmpz_t den, slong len,
                                const fmpz_t c)
{
    fmpz_t gcd;  /* GCD( den, c ) */

    if (fmpz_is_zero(c))
    {
        _fmpz_vec_zero(rpoly, len);
        fmpz_one(rden);
        return;
    }

    fmpz_init(gcd);
    fmpz_one(gcd);
    if (!fmpz_is_one(c))
        fmpz_gcd(gcd, c, den);
    if (fmpz_is_one(gcd))
    {
        _fmpz_vec_scalar_mul_fmpz(rpoly, poly, len, c);
        fmpz_set(rden, den);
    }
    else
    {
        fmpz_t c2;
        fmpz_init(c2);
        fmpz_divexact(c2, c, gcd);
        _fmpz_vec_scalar_mul_fmpz(rpoly, poly, len, c2);
        fmpz_divexact(rden, den, gcd);
        fmpz_clear(c2);
    }
    fmpz_clear(gcd);
}

void fmpq_poly_scalar_mul_fmpz(fmpq_poly_t rop, const fmpq_poly_t op, const fmpz_t c)
{
    if (fmpz_is_zero(c) || fmpq_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
        return;
    }
    
    fmpq_poly_fit_length(rop, op->length);
    _fmpq_poly_set_length(rop, op->length);
    
    _fmpq_poly_scalar_mul_fmpz(rop->coeffs, rop->den, 
                               op->coeffs, op->den, op->length, c);
}

