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

void _fmpz_poly_rem_powers_precomp(fmpz * A, slong m, 
                                 const fmpz * B, slong n, fmpz ** const powers)
{
   slong i;
   fmpz * prod;

   if (m >= 2*n)
   {
      _fmpz_poly_rem_powers_precomp(A + n, m - n, B, n, powers);
      m = 2*n - 1;
      while (m && fmpz_is_zero(A + m - 1)) m--;
   }

   if (m < n)
      return;

   prod = _fmpz_vec_init(n - 1);
             
   for (i = n - 1; i < m; i++)
   {
      _fmpz_vec_scalar_mul_fmpz(prod, powers[i], n - 1, A + i);
      _fmpz_poly_add(A, A, n - 1, prod, n - 1);
   }

   _fmpz_vec_clear(prod, n - 1);
}

void 
fmpz_poly_rem_powers_precomp(fmpz_poly_t R, const fmpz_poly_t A, 
                    const fmpz_poly_t B, const fmpz_poly_powers_precomp_t B_inv)
{
    fmpz_poly_t tR;
    fmpz *r;
    slong len1 = A->length, len2 = B->length;
    
    if (len1 < len2)
    {
        fmpz_poly_set(R, A);
        return;
    }

    if (R == B)
    {
        fmpz_poly_init2(tR, len1);
        r = tR->coeffs;
    }
    else
    {
        fmpz_poly_fit_length(R, len1);
        r = R->coeffs;
    }

    if (R == B || R != A)
       _fmpz_vec_set(r, A->coeffs, len1);

    _fmpz_poly_rem_powers_precomp(r, len1, B->coeffs, len2, B_inv->powers);

    if (R == B)
    {
        _fmpz_poly_set_length(tR, len2 - 1);
        fmpz_poly_swap(tR, R);
        fmpz_poly_clear(tR);
    }
    else
        _fmpz_poly_set_length(R, len2 - 1);

    _fmpz_poly_normalise(R);
}
