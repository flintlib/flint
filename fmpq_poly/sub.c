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

void _fmpq_poly_sub_can(fmpz * rpoly, fmpz_t rden, 
                    const fmpz * poly1, const fmpz_t den1, slong len1, 
                    const fmpz * poly2, const fmpz_t den2, slong len2, int can)
{
    slong max = FLINT_MAX(len1, len2);
    slong min = FLINT_MIN(len1, len2);
    
    fmpz_t d;

    if (fmpz_equal(den1, den2))
    {
        _fmpz_poly_sub(rpoly, poly1, len1, poly2, len2);

        if (fmpz_is_one(den1) || !can)
            fmpz_set(rden, den1);
        else if (can)
        {
            fmpz_init(d);
            _fmpz_vec_content_chained(d, rpoly, max, den1);

            if (fmpz_is_one(d))
                  fmpz_set(rden, den1);
            else
            {
                _fmpz_vec_scalar_divexact_fmpz(rpoly, rpoly, max, d);
                fmpz_divexact(rden, den1, d);
            }

            fmpz_clear(d);
        }

        return;
    }

    fmpz_init(d);
    fmpz_one(d);
    if (!fmpz_is_one(den1) && !fmpz_is_one(den2))
        fmpz_gcd(d, den1, den2);
    
    if (fmpz_is_one(d))
    {
        _fmpz_vec_scalar_mul_fmpz(rpoly, poly1, len1, den2);
        _fmpz_vec_scalar_submul_fmpz(rpoly, poly2, min, den1);
        if (len1 < len2)
        {
            _fmpz_vec_scalar_mul_fmpz(rpoly + min, poly2 + min, max - min, den1);
            _fmpz_vec_neg(rpoly + min, rpoly + min, max - min);
        }
        fmpz_mul(rden, den1, den2);
    }
    else
    {
        fmpz_t den11;
        fmpz_t den22;
        fmpz_init(den11);
        fmpz_init(den22);
        fmpz_divexact(den11, den1, d);
        fmpz_divexact(den22, den2, d);
        
        _fmpz_vec_scalar_mul_fmpz(rpoly, poly1, len1, den22);
        _fmpz_vec_scalar_submul_fmpz(rpoly, poly2, len2, den11);
        if (len1 < len2)
        {
            _fmpz_vec_scalar_mul_fmpz(rpoly + min, poly2 + min, max - min, den11);
            _fmpz_vec_neg(rpoly + min, rpoly + min, max - min);
        }
        
        if (_fmpz_vec_is_zero(rpoly, max))
            fmpz_one(rden);
        else
        {
            if (can)
            {
               fmpz_t e;
               fmpz_init(e);
               _fmpz_vec_content_chained(e, rpoly, max, d);
            
               if (fmpz_is_one(e))
                  fmpz_mul(rden, den1, den22);
               else
               {
                  _fmpz_vec_scalar_divexact_fmpz(rpoly, rpoly, max, e);
                  fmpz_divexact(den11, den1, e);
                  fmpz_mul(rden, den11, den22);
               }
               fmpz_clear(e);
            } else
               fmpz_mul(rden, den1, den22);
        }
        fmpz_clear(den11);
        fmpz_clear(den22);
    }
    fmpz_clear(d);
}

void _fmpq_poly_sub(fmpz * rpoly, fmpz_t rden, 
                    const fmpz * poly1, const fmpz_t den1, slong len1, 
                    const fmpz * poly2, const fmpz_t den2, slong len2)
{
   _fmpq_poly_sub_can(rpoly, rden, poly1, den1, len1, poly2, den2, len2, 1);
}

void fmpq_poly_sub_can(fmpq_poly_t res, const fmpq_poly_t poly1, 
                                        const fmpq_poly_t poly2, int can)
{
    slong len1, len2, max;

    if (poly1 == poly2)
    {
        fmpq_poly_zero(res);
        return;
    }
    
    len1 = poly1->length;
    len2 = poly2->length;
    max  = FLINT_MAX(poly1->length, poly2->length);
    fmpq_poly_fit_length(res, max);
    
    if (res != poly2)
        _fmpq_poly_sub_can(res->coeffs, res->den, 
                       poly1->coeffs, poly1->den, len1, 
                       poly2->coeffs, poly2->den, len2, can);
    else
    {
        _fmpq_poly_sub_can(res->coeffs, res->den, 
                       poly2->coeffs, poly2->den, len2, 
                       poly1->coeffs, poly1->den, len1, can);
        _fmpz_vec_neg(res->coeffs, res->coeffs, max);
    }
    
    _fmpq_poly_set_length(res, max);
    _fmpq_poly_normalise(res);
}

void fmpq_poly_sub(fmpq_poly_t res, const fmpq_poly_t poly1, 
                                        const fmpq_poly_t poly2)
{
    fmpq_poly_sub_can(res, poly1, poly2, 1);
}

void fmpq_poly_sub_si(fmpq_poly_t res, const fmpq_poly_t poly, slong c)
{
    if (c == 0)
    {
        fmpq_poly_set(res, poly);
    }
    else if (poly->length == 0)
    {
        fmpq_poly_set_si(res, c);
        fmpz_neg(res->coeffs + 0, res->coeffs + 0);
    }
    else
    {
        fmpz_t p, q;
        fmpz_init_set_si(p, c);
        *q = 1;
        fmpq_poly_fit_length(res, poly->length);
        _fmpq_poly_set_length(res, poly->length);
        _fmpq_poly_sub(res->coeffs, res->den, poly->coeffs, poly->den, poly->length, p, q, 1);
        _fmpq_poly_normalise(res);
        fmpz_clear(p);
    }
}

void fmpq_poly_si_sub(fmpq_poly_t res, slong c, const fmpq_poly_t poly)
{
    if (c == 0)
    {
        fmpq_poly_neg(res, poly);
    }
    else if (poly->length == 0)
    {
        fmpq_poly_set_si(res, c);
    }
    else
    {
        fmpz_t p, q;
        fmpz_init_set_si(p, c);
        *q = 1;
        fmpq_poly_fit_length(res, poly->length);
        _fmpq_poly_set_length(res, poly->length);
        _fmpq_poly_sub(res->coeffs, res->den, poly->coeffs, poly->den, poly->length, p, q, 1);
        _fmpq_poly_normalise(res);
        _fmpz_vec_neg(res->coeffs, res->coeffs, res->length);
        fmpz_clear(p);
    }
}

void fmpq_poly_sub_fmpz(fmpq_poly_t res, const fmpq_poly_t poly, const fmpz_t c)
{
    if (fmpz_is_zero(c))
    {
        fmpq_poly_set(res, poly);
    }
    else if (poly->length == 0)
    {
        fmpq_poly_set_fmpz(res, c);
        fmpz_neg(res->coeffs + 0, res->coeffs + 0);
    }
    else
    {
        fmpz_t q;
        *q = 1;
        fmpq_poly_fit_length(res, poly->length);
        _fmpq_poly_set_length(res, poly->length);
        _fmpq_poly_sub(res->coeffs, res->den, poly->coeffs, poly->den, poly->length, c, q, 1);
        _fmpq_poly_normalise(res);
    }
}

void fmpq_poly_fmpz_sub(fmpq_poly_t res, const fmpz_t c, const fmpq_poly_t poly)
{
    if (fmpz_is_zero(c))
    {
        fmpq_poly_neg(res, poly);
    }
    else if (poly->length == 0)
    {
        fmpq_poly_set_fmpz(res, c);
    }
    else
    {
        fmpz_t q;
        *q = 1;
        fmpq_poly_fit_length(res, poly->length);
        _fmpq_poly_set_length(res, poly->length);
        _fmpq_poly_sub(res->coeffs, res->den, poly->coeffs, poly->den, poly->length, c, q, 1);
        _fmpq_poly_normalise(res);
        _fmpz_vec_neg(res->coeffs, res->coeffs, res->length);
    }
}

void fmpq_poly_sub_fmpq(fmpq_poly_t res, const fmpq_poly_t poly, const fmpq_t c)
{
    if (fmpq_is_zero(c))
    {
        fmpq_poly_set(res, poly);
    }
    else if (poly->length == 0)
    {
        fmpq_poly_set_fmpq(res, c);
        fmpz_neg(res->coeffs + 0, res->coeffs + 0);
    }
    else
    {
        fmpq_poly_fit_length(res, poly->length);
        _fmpq_poly_set_length(res, poly->length);
        _fmpq_poly_sub(res->coeffs, res->den, poly->coeffs, poly->den, poly->length, fmpq_numref(c), fmpq_denref(c), 1);
        _fmpq_poly_normalise(res);
    }
}

void fmpq_poly_fmpq_sub(fmpq_poly_t res, const fmpq_t c, const fmpq_poly_t poly)
{
    if (fmpq_is_zero(c))
    {
        fmpq_poly_neg(res, poly);
    }
    else if (poly->length == 0)
    {
        fmpq_poly_set_fmpq(res, c);
    }
    else
    {
        fmpq_poly_fit_length(res, poly->length);
        _fmpq_poly_set_length(res, poly->length);
        _fmpq_poly_sub(res->coeffs, res->den, poly->coeffs, poly->den, poly->length, fmpq_numref(c), fmpq_denref(c), 1);
        _fmpq_poly_normalise(res);
        _fmpz_vec_neg(res->coeffs, res->coeffs, res->length);
    }
}
