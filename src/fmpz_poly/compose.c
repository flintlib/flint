/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2016 Shivin Srivastava

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

void
_fmpz_poly_compose(fmpz * res, const fmpz * poly1, slong len1, 
                               const fmpz * poly2, slong len2)
{
    if (len1 == 1)
        fmpz_set(res, poly1);
    else if (len2 == 1)
        _fmpz_poly_evaluate_fmpz(res, poly1, len1, poly2);
    else if (len1 <= 4)
        _fmpz_poly_compose_horner(res, poly1, len1, poly2, len2);
    else if (len2 == 2)
    {
        slong i;
        _fmpz_vec_set(res, poly1, len1);
        _fmpz_poly_taylor_shift(res, poly2, len1);
        
        if (fmpz_equal_si(poly2 + 1, -1))
        {
            for (i = 1; i < len1; i += 2)
                fmpz_neg(res + i, res + i);
            return;
        }
        else if (!fmpz_is_one(poly2 + 1))
        {
            fmpz_t temp;
            fmpz_init(temp);
            fmpz_one(temp);
        
            for (i = 0; i < len1; i++)
            {
                fmpz_mul(res + i, res + i, temp);
                /* no need to reverse signs manually as if poly2 + 1 negative 
                then signs alternate automatically*/
                fmpz_mul(temp, temp, poly2 + 1);
            }
            fmpz_clear(temp);
        }
        return; 
    }
    else
        _fmpz_poly_compose_divconquer(res, poly1, len1, poly2, len2);
}

void
fmpz_poly_compose(fmpz_poly_t res, 
                  const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    const slong len1 = poly1->length;
    const slong len2 = poly2->length;
    slong lenr;
    
    if (len1 == 0)
    {
        fmpz_poly_zero(res);
        return;
    }
    if (len1 == 1 || len2 == 0)
    {
        fmpz_poly_set_fmpz(res, poly1->coeffs);
        return;
    }
  
    lenr = (len1 - 1) * (len2 - 1) + 1;
    
    if (res != poly1 && res != poly2)
    {
        fmpz_poly_fit_length(res, lenr);
        _fmpz_poly_compose(res->coeffs, poly1->coeffs, len1, 
                                        poly2->coeffs, len2);
        _fmpz_poly_set_length(res, lenr);
        _fmpz_poly_normalise(res);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, lenr);
        _fmpz_poly_compose(t->coeffs, poly1->coeffs, len1,
                                      poly2->coeffs, len2);
        _fmpz_poly_set_length(t, lenr);
        _fmpz_poly_normalise(t);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
}
