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

void
_fmpq_poly_compose(fmpz * res, fmpz_t den, const fmpz * poly1, const fmpz_t den1, 
                   slong len1, const fmpz * poly2, const fmpz_t den2, slong len2)
{
    if (*den2 == WORD(1))
    {
        _fmpz_poly_compose(res, poly1, len1, poly2, len2);
        fmpz_set(den, den1);
        _fmpq_poly_canonicalise(res, den, (len1 - WORD(1)) * (len2 - WORD(1)) + WORD(1));
    }
    else
    {
        fmpz_t one;
        fmpz * v = _fmpz_vec_init(len1);
        fmpz_init(one);
        fmpz_one(one);
        
        _fmpq_poly_rescale(v, den, poly1, den1, len1, one, den2);
        _fmpz_poly_compose(res, v, len1, poly2, len2);
        _fmpq_poly_canonicalise(res, den, (len1 - WORD(1)) * (len2 - WORD(1)) + WORD(1));
        
        fmpz_clear(one);
        _fmpz_vec_clear(v, len1);
    }
}

void
fmpq_poly_compose(fmpq_poly_t res, 
                              const fmpq_poly_t poly1, const fmpq_poly_t poly2)
{
    const slong len1 = poly1->length;
    const slong len2 = poly2->length;
    slong lenr;
    
    if (len1 == WORD(0))
    {
        fmpq_poly_zero(res);
        return;
    }
    if (len1 == WORD(1) || len2 == WORD(0))
    {
        fmpq_poly_fit_length(res, 1);
        fmpz_set(res->coeffs, poly1->coeffs);
        fmpz_set(res->den, poly1->den);
        {
            fmpz_t d;
            fmpz_init(d);
            fmpz_gcd(d, res->coeffs, res->den);
            if (*d != WORD(1))
            {
                fmpz_divexact(res->coeffs, res->coeffs, d);
                fmpz_divexact(res->den, res->den, d);
            }
            fmpz_clear(d);
        }
        _fmpq_poly_set_length(res, 1);
        _fmpq_poly_normalise(res);
        return;
    }
    
    lenr = (len1 - WORD(1)) * (len2 - WORD(1)) + WORD(1);
    
    if ((res != poly1) && (res != poly2))
    {
        fmpq_poly_fit_length(res, lenr);
        _fmpq_poly_compose(res->coeffs, res->den, 
                           poly1->coeffs, poly1->den, len1, 
                           poly2->coeffs, poly2->den, len2);
        _fmpq_poly_set_length(res, lenr);
        _fmpq_poly_normalise(res);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, lenr);
        _fmpq_poly_compose(t->coeffs, t->den, 
                           poly1->coeffs, poly1->den, len1,
                           poly2->coeffs, poly2->den, len2);
        _fmpq_poly_set_length(t, lenr);
        _fmpq_poly_normalise(t);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
    }
}
