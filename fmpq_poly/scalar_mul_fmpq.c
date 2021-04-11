/*
    Copyright (C) 2011 Sebastian Pancratz

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

void _fmpq_poly_scalar_mul_fmpq(fmpz * rpoly, fmpz_t rden, 
                                const fmpz * poly, const fmpz_t den, slong len, 
                                const fmpz_t r, const fmpz_t s)
{
    fmpz_t gcd1;  /* GCD( poly, s ) */
    fmpz_t gcd2;  /* GCD( r, den )  */

    if (fmpz_is_zero(r))
    {
        _fmpz_vec_zero(rpoly, len);
        fmpz_one(rden);
        return;
    }

    fmpz_init(gcd1);
    fmpz_init(gcd2);
    fmpz_one(gcd1);
    fmpz_one(gcd2);
    if (!fmpz_is_one(s))
        _fmpz_vec_content_chained(gcd1, poly, len, s);
    if (!fmpz_is_one(den) && !fmpz_is_one(r))
        fmpz_gcd(gcd2, r, den);
    
    if (fmpz_is_one(gcd1))
    {
        if (fmpz_is_one(gcd2))
        {
            _fmpz_vec_scalar_mul_fmpz(rpoly, poly, len, r);
            fmpz_mul(rden, den, s);
        }
        else
        {
            fmpz_t r2;
            fmpz_init(r2);
            fmpz_divexact(r2, r, gcd2);
            _fmpz_vec_scalar_mul_fmpz(rpoly, poly, len, r2);
            fmpz_divexact(rden, den, gcd2);
            fmpz_mul(rden, rden, s);
            fmpz_clear(r2);
        }
    }
    else
    {
        fmpz_t s2;
        fmpz_init(s2);
        fmpz_divexact(s2, s, gcd1);
        if (fmpz_is_one(gcd2))
        {
            _fmpz_vec_scalar_divexact_fmpz(rpoly, poly, len, gcd1);
            _fmpz_vec_scalar_mul_fmpz(rpoly, rpoly, len, r);
            fmpz_mul(rden, den, s2);
        }
        else
        {
            fmpz_t r2;
            fmpz_init(r2);
            fmpz_divexact(r2, r, gcd2);
            _fmpz_vec_scalar_divexact_fmpz(rpoly, poly, len, gcd1);
            _fmpz_vec_scalar_mul_fmpz(rpoly, rpoly, len, r2);
            fmpz_divexact(rden, den, gcd2);
            fmpz_mul(rden, rden, s2);
            fmpz_clear(r2);
        }
        fmpz_clear(s2);
    }
    
    fmpz_clear(gcd1);
    fmpz_clear(gcd2);
}

void fmpq_poly_scalar_mul_fmpq(fmpq_poly_t rop, const fmpq_poly_t op, const fmpq_t c)
{
    if (fmpz_is_one(fmpq_denref(c)))
    {
        fmpq_poly_scalar_mul_fmpz(rop, op, fmpq_numref(c));
    }
    else if (fmpq_is_zero(c) || fmpq_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
    }
    else
    {
        fmpq_poly_fit_length(rop, op->length);
        _fmpq_poly_set_length(rop, op->length);
        
        _fmpq_poly_scalar_mul_fmpq(rop->coeffs, rop->den, 
                                   op->coeffs, op->den, op->length, 
                                   fmpq_numref(c), fmpq_denref(c));
    }
}

