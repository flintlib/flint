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
#include "fmpq.h"
#include "fmpq_poly.h"

void _fmpq_poly_scalar_div_fmpq(fmpz * rpoly, fmpz_t rden,
                                const fmpz * poly, const fmpz_t den, slong len,
                                const fmpz_t r, const fmpz_t s)
{
    fmpz_t gcd1;  /* GCD( poly, r ) */
    fmpz_t gcd2;  /* GCD( s, den )  */
    fmpz_init(gcd1);
    fmpz_init(gcd2);
    fmpz_one(gcd1);
    fmpz_one(gcd2);
    if (!fmpz_is_one(r))
        _fmpz_vec_content_chained(gcd1, poly, len, r);
    if (!fmpz_is_one(den) && !fmpz_is_one(s))
        fmpz_gcd(gcd2, s, den);

    if (fmpz_is_one(gcd1))
    {
        if (fmpz_is_one(gcd2))
        {
            _fmpz_vec_scalar_mul_fmpz(rpoly, poly, len, s);
            fmpz_mul(rden, den, r);
        }
        else
        {
            fmpz_t s2;
            fmpz_init(s2);
            fmpz_divexact(s2, s, gcd2);
            _fmpz_vec_scalar_mul_fmpz(rpoly, poly, len, s2);
            fmpz_divexact(rden, den, gcd2);
            fmpz_mul(rden, rden, r);
            fmpz_clear(s2);
        }
    }
    else
    {
        fmpz_t r2;
        fmpz_init(r2);
        fmpz_divexact(r2, r, gcd1);
        if (fmpz_is_one(gcd2))
        {
            _fmpz_vec_scalar_divexact_fmpz(rpoly, poly, len, gcd1);
            _fmpz_vec_scalar_mul_fmpz(rpoly, rpoly, len, s);
            fmpz_mul(rden, den, r2);
        }
        else
        {
            fmpz_t s2;
            fmpz_init(s2);
            fmpz_divexact(s2, s, gcd2);
            _fmpz_vec_scalar_divexact_fmpz(rpoly, poly, len, gcd1);
            _fmpz_vec_scalar_mul_fmpz(rpoly, rpoly, len, s2);
            fmpz_divexact(rden, den, gcd2);
            fmpz_mul(rden, rden, r2);
            fmpz_clear(s2);
        }
        fmpz_clear(r2);
    }

    if (_fmpz_vec_is_zero(rpoly, len))
        fmpz_one(rden);
    if (fmpz_sgn(rden) < 0)
    {
        _fmpz_vec_neg(rpoly, rpoly, len);
        fmpz_neg(rden, rden);
    }

    fmpz_clear(gcd1);
    fmpz_clear(gcd2);
}

void fmpq_poly_scalar_div_fmpq(fmpq_poly_t rop, const fmpq_poly_t op, const fmpq_t c)
{
    if (fmpq_is_zero(c))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_poly_scalar_div_fmpq). Division by zero.\n");
    }

    if (fmpq_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
    }
    else
    {
        fmpq_poly_fit_length(rop, op->length);
        _fmpq_poly_set_length(rop, op->length);

        _fmpq_poly_scalar_div_fmpq(rop->coeffs, rop->den,
                                   op->coeffs, op->den, op->length,
                                   fmpq_numref(c), fmpq_denref(c));
    }
}
