/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2014 William Hart

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

int _fmpq_poly_equal_trunc(const fmpz * poly1, const fmpz_t den1, slong len1, 
                           const fmpz * poly2, const fmpz_t den2, slong len2, slong n)
{
    int res = 1;
    slong i;

    if (n < 0)
       n = 0;

    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);

    if (len1 > len2) 
    {
        for (i = len2; i < len1; i++)
        {
           if (!fmpz_is_zero(poly1 + i))
              return 0;
        }
        
        len1 = len2;
    } else if (len2 > len1)
    {
        for (i = len1; i < len2; i++)
        {
           if (!fmpz_is_zero(poly2 + i))
              return 0;
        }
    }

    if (fmpz_equal(den1, den2))
        return (_fmpz_vec_equal(poly1, poly2, len1));
    else
    {
       fmpz_t p1, p2, d, d1, d2;

       fmpz_init(d);
       fmpz_init(p1);
       fmpz_init(p2);
       fmpz_init(d1);
       fmpz_init(d2);

       fmpz_gcd(d, den1, den2);

       if (!fmpz_is_one(d))
       {
          fmpz_divexact(d1, den1, d);
          fmpz_divexact(d2, den2, d);
       } else
       {
          fmpz_set(d1, den1);
          fmpz_set(d2, den2);
       }

       for (i = 0; i < len1; i++)
       {
          fmpz_mul(p1, poly1 + i, d2);
          fmpz_mul(p2, poly2 + i, d1);
          if (!fmpz_equal(p1, p2))
          {
             res = 0;
             break;
          }
       }
       
       fmpz_clear(d1);
       fmpz_clear(d2);
       fmpz_clear(p1);
       fmpz_clear(p2);
       fmpz_clear(d);
    }

    return res;
}

int fmpq_poly_equal_trunc(const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n)
{
   return _fmpq_poly_equal_trunc(poly1->coeffs, poly1->den, poly1->length,
                                 poly2->coeffs, poly2->den, poly2->length, n);
}
