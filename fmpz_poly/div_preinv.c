/*
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void
_fmpz_poly_div_preinv(fmpz * Q, const fmpz * A, slong len1_in, 
                                const fmpz * B, const fmpz * B_inv, slong len2)
{
   slong len1 = len1_in;
   slong n = len1 - len2 + 1;
   fmpz * A_rev;
   fmpz * a;
   
   if (n > len2)
   {
      a = _fmpz_vec_init(len1_in);
      _fmpz_vec_set(a, A, len1_in);
      
      do {
         slong start = n - len2;
         _fmpz_poly_divrem_preinv(Q + start, a + start, len1 - start, 
                                                               B, B_inv, len2);
         n -= len2;
         len1 -= len2;
      } while (n > len2);
   } else
      a = (fmpz *) A;

   A_rev = _fmpz_vec_init(len1);
   
   _fmpz_poly_reverse(A_rev, a, len1, len1);
   _fmpz_poly_mullow(Q, A_rev, len1, B_inv, len2, n);
   _fmpz_poly_reverse(Q, Q, n, n);

   if (a != A)
      _fmpz_vec_clear(a, len1_in);
   _fmpz_vec_clear(A_rev, len1);
}

void
fmpz_poly_div_preinv(fmpz_poly_t Q, const fmpz_poly_t A, 
                                  const fmpz_poly_t B, const fmpz_poly_t B_inv)
{
    fmpz_poly_t tQ;
    fmpz *q;
    slong len1 = A->length, len2 = B_inv->length;
    slong qlen = len1 - len2 + 1;

    if (len1 < len2)
    {
        fmpz_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B || Q == B_inv)
    {
        fmpz_poly_init2(tQ, qlen);
        q = tQ->coeffs;
    }
    else
    {
        fmpz_poly_fit_length(Q, qlen);
        q = Q->coeffs;
    }
    
    _fmpz_poly_div_preinv(q, A->coeffs, len1, B->coeffs, B_inv->coeffs, len2);

    if (Q == A || Q == B || Q == B_inv)
    {
        _fmpz_poly_set_length(tQ, qlen);
        fmpz_poly_swap(tQ, Q);
        fmpz_poly_clear(tQ);
    }
    else
        _fmpz_poly_set_length(Q, qlen);

    /* no need to normalise */
}
