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
_fmpq_poly_rescale(fmpz * res, fmpz_t denr, 
                   const fmpz * poly, const fmpz_t den, slong len, 
                   const fmpz_t xnum, const fmpz_t xden)
{
    if (len < WORD(2))
    {
        if (res != poly)
        {
            _fmpz_vec_set(res, poly, len);
            fmpz_set(denr, den);
        }
    }
    else
    {
        slong i;
        fmpz_t t;
        
        fmpz_init(t);
        
        fmpz_one(t);
        fmpz_set(res, poly);
        for (i = WORD(1); i < len; i++)
        {
            fmpz_mul(t, t, xnum);
            fmpz_mul(res + i, poly + i, t);
        }
        fmpz_one(t);
        for (i = len - WORD(2); i >= WORD(0); i--)
        {
            fmpz_mul(t, t, xden);
            fmpz_mul(res + i, res + i, t);
        }
        fmpz_mul(denr, den, t);
        
        fmpz_clear(t);
        
        _fmpq_poly_canonicalise(res, denr, len);
    }
}

void fmpq_poly_rescale(fmpq_poly_t res, const fmpq_poly_t poly, const fmpq_t x)
{
    if (fmpq_is_zero(x))
    {
        fmpq_poly_zero(res);
    }
    else if (poly->length < WORD(2))
    {
        fmpq_poly_set(res, poly);
    }
    else
    {
        fmpq_poly_fit_length(res, poly->length);
        _fmpq_poly_rescale(res->coeffs, res->den, 
                           poly->coeffs, poly->den, poly->length, 
                           fmpq_numref(x), fmpq_denref(x));
        _fmpq_poly_set_length(res, poly->length);
    }
}
