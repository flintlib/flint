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
#include "fmpz_poly.h"
#include "fmpq_poly.h"

void _fmpq_poly_mul(fmpz * rpoly, fmpz_t rden, 
                    const fmpz * poly1, const fmpz_t den1, slong len1, 
                    const fmpz * poly2, const fmpz_t den2, slong len2)
{
    fmpz_t gcd1;  /* GCD( poly1, den2 ) */
    fmpz_t gcd2;  /* GCD( poly2, den1 ) */

    if (poly1 == poly2 && len1 == len2)
    {
        _fmpz_poly_sqr(rpoly, poly1, len1);
        fmpz_mul(rden, den1, den2);
        return;
    }

    fmpz_init(gcd1);
    fmpz_init(gcd2);
    fmpz_one(gcd1);
    fmpz_one(gcd2);
    
    if (!fmpz_is_one(den2))
        _fmpz_vec_content_chained(gcd1, poly1, len1, den2);

    if (!fmpz_is_one(den1))
        _fmpz_vec_content_chained(gcd2, poly2, len2, den1);

    /*
       TODO:  If gcd1 and gcd2 are very large compared to the degrees of 
       poly1 and poly2, we might want to create copies of the polynomials 
       and divide out the common factors *before* the multiplication.
     */
    
    _fmpz_poly_mul(rpoly, poly1, len1, poly2, len2);
    fmpz_mul(rden, den1, den2);

    if (!fmpz_is_one(gcd1) || !fmpz_is_one(gcd2))
    {
        fmpz_t g;
        fmpz_init(g);
        fmpz_mul(g, gcd1, gcd2);
        _fmpz_vec_scalar_divexact_fmpz(rpoly, rpoly, len1 + len2 - 1, g);
        fmpz_divexact(rden, rden, g);
        fmpz_clear(g);
    }
    
    fmpz_clear(gcd1);
    fmpz_clear(gcd2);
}

void fmpq_poly_mul(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2)
{
    slong len;
    
    if (poly1->length == 0 || poly2->length == 0)
    {
        fmpq_poly_zero(res);
        return;
    }

    len = poly1->length + poly2->length - 1;

    if (res == poly2 || res == poly1)
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, len);
        fmpq_poly_mul(t, poly1, poly2);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
        return;
    }
    
    fmpq_poly_fit_length(res, len);

    if (poly1->length >= poly2->length)
        _fmpq_poly_mul(res->coeffs, res->den, 
                       poly1->coeffs, poly1->den, poly1->length,
                       poly2->coeffs, poly2->den, poly2->length);
    else
        _fmpq_poly_mul(res->coeffs, res->den, 
                       poly2->coeffs, poly2->den, poly2->length,
                       poly1->coeffs, poly1->den, poly1->length);

    _fmpq_poly_set_length(res, len);
}

