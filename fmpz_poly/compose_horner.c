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

void
_fmpz_poly_compose_horner(fmpz * res, const fmpz * poly1, slong len1, 
                                      const fmpz * poly2, slong len2)
{
    if (len1 == 1)
    {
        fmpz_set(res, poly1);
    }
    else
    {
        const slong alloc = (len1 - 1) * (len2 - 1) + 1;

        slong i = len1 - 1, lenr;
        fmpz * t = _fmpz_vec_init(alloc);
        
        /*
           Perform the first two steps as one, 
             "res = a(m) * poly2 + a(m-1)".
         */
        {
            lenr = len2;
            _fmpz_vec_scalar_mul_fmpz(res, poly2, len2, poly1 + i);
            i--;
            fmpz_add(res, res, poly1 + i);
        }
        while (i--)
        {
            _fmpz_poly_mul(t, res, lenr, poly2, len2);
            lenr += len2 - 1;
            _fmpz_poly_add(res, t, lenr, poly1 + i, 1);
        }
        _fmpz_vec_clear(t, alloc);
    }
}

void
fmpz_poly_compose_horner(fmpz_poly_t res, 
                         const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    const slong len1 = poly1->length;
    const slong len2 = poly2->length;
    
    if (len1 == 0)
    {
        fmpz_poly_zero(res);
    }
    else if (len1 == 1 || len2 == 0)
    {
        fmpz_poly_set_fmpz(res, poly1->coeffs);
    }
    else
    {
        const slong lenr = (len1 - 1) * (len2 - 1) + 1;
        
        if ((res != poly1) && (res != poly2))
        {
            fmpz_poly_fit_length(res, lenr);
            _fmpz_poly_compose_horner(res->coeffs, poly1->coeffs, len1, 
                                                   poly2->coeffs, len2);
            _fmpz_poly_set_length(res, lenr);
        }
        else
        {
            fmpz_poly_t t;
            fmpz_poly_init2(t, lenr);
            _fmpz_poly_compose_horner(t->coeffs, poly1->coeffs, len1,
                                                 poly2->coeffs, len2);
            _fmpz_poly_set_length(t, lenr);
            fmpz_poly_swap(res, t);
            fmpz_poly_clear(t);
        }
        _fmpz_poly_normalise(res);
    }
}
